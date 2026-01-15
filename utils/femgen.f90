!===============================================================================
! FEMGEN - Finite Element Modal Generator for MBDyn
!===============================================================================
!
! PURPOSE
!   Converts modal analysis results from MSC Nastran SOL 103 (normal modes)
!   into a modal-data file (.fem) suitable for MBDyn's modal joint element.
!
! INPUT FILES (must be in current directory)
!   mbdyn.mat : Nastran OUTPUT4 binary file containing matrices:
!               - MHH    : Modal mass matrix (n_modes x n_modes)
!               - KHH    : Modal stiffness matrix (n_modes x n_modes)
!               - PHIX   : Mode shape matrix (6*n_nodes x n_modes)
!               - LUMPMS : Lumped mass diagonal (6*n_nodes x 1)
!
!   mbdyn.tab : Nastran OUTPUT2 binary file containing tables:
!               - GPL    : Grid Point List (node IDs)
!               - BGPDT  : Basic Grid Point Definition Table (coordinates)
!
! OUTPUT FILE
!   <name>.fem : MBDyn modal data file with record groups:
!                1. Header (revision, counts)
!                2. Node list
!                3. Initial modal displacements (optional)
!                4. Initial modal velocities (optional)
!                5-7. Nodal X, Y, Z coordinates
!                8. Mode shapes
!                9. Modal mass matrix
!                10. Modal stiffness matrix
!                11. Lumped mass diagonal
!
! PUBLIC INTERFACE
!   integer(c_int) function femgen(outname, 
!                                    is_modal_displacement, 
!                                    is_modal_velocity, 
!                                    mass_direction) bind(C)
!     outname : Output filename (C string, max 72 chars)
!     is_modal_displacement    : Include initial modal displacements (0=no, 1=yes)
!     is_modal_velocity        : Include initial modal velocities (0=no, 1=yes)
!     mass_direction           : Lumped mass treatment (-1=auto, 0=none, 1-3=use component)
!     Returns                  : 0 on success, non-zero on failure
!
! AUTHORS
!   Giuseppe Quaranta: <giuseppe.quaranta@polimi.it>
!   Alessandro Cocco : <alessandro.cocco@polimi.it> 
!
!===============================================================================

module femgen_mod
  use, intrinsic :: iso_fortran_env, only: &
    int32, real32, real64
  use, intrinsic :: iso_c_binding, only: &
    c_int, c_int32_t, c_char, c_null_char, c_float, c_double
  
  implicit none
  private

  public :: femgen_c_wrapper

  !---------------------------------------------------------------------------
  ! Module Constants
  !---------------------------------------------------------------------------
  !> Tolerance for snapping small values to zero in output
  real(real64), parameter :: FEMGEN_ZERO_TOL = 1.0e-16_real64
  
  !> Maximum number of FEM nodes supported
  integer, parameter :: MAX_NODES = 20000
  
  !> Maximum words per OUTPUT4 matrix record
  integer, parameter :: OUTPUT4_MAX_RECORD_WORDS = 2000000
  
  !> Maximum words per OUTPUT2 data block
  integer, parameter :: OUTPUT2_MAX_BLOCK_WORDS = 5000000

  !> Maximum length of output filename (including null terminator) 
  integer, parameter :: MAX_OUTNAME_LENGTH = 73
  
  !> Default input filenames (Nastran convention)
  character(len=9), parameter :: FILE_MBDYN_MAT = 'mbdyn.mat' 
  character(len=9), parameter :: FILE_MBDYN_TAB = 'mbdyn.tab'

  !---------------------------------------------------------------------------
  ! Derived Types
  !---------------------------------------------------------------------------
  
  !> Container for complete modal model data
  !>
  !> Holds all information extracted from Nastran files needed to generate
  !> the MBDyn .fem output file.
  type :: femgen_model
    integer :: n_modes = 0                       !< Number of vibration modes
    integer :: n_nodes = 0                       !< Number of FEM nodes
    integer :: phix_nrow = 0                     !< PHIX (eigenvectors) row dimension
    integer :: phix_ncol = 0                     !< PHIX (eigenvectors) column dimension
    integer(int32), allocatable :: node_ids(:)   !< External node ID labels [n_nodes]
    integer(int32), allocatable :: internal_ids(:) !< Internal NASTRAN IDs [n_nodes]
    real(real64), allocatable :: xg(:)           !< Global X coordinates [n_nodes]
    real(real64), allocatable :: yg(:)           !< Global Y coordinates [n_nodes]
    real(real64), allocatable :: zg(:)           !< Global Z coordinates [n_nodes]
    real(real64), allocatable :: mhh(:,:)        !< Modal mass matrix [n_modes, n_modes]
    real(real64), allocatable :: khh(:,:)        !< Modal stiffness matrix [n_modes, n_modes]
    real(real64), allocatable :: lump_vec(:)     !< Lumped mass diagonal [6*n_nodes]
    real(real64), allocatable :: phix(:,:)       !< Mode shape matrix [6*n_nodes, n_modes] 
  end type femgen_model

  !> Nastran OUTPUT2 binary file reader
  !>
  !> Provides buffered reading of Nastran OUTPUT2 format files, handling
  !> the complex record structure with variable-length data blocks.
  type :: output2_reader
    integer :: file_unit = -1                 !< Fortran unit number for input file
    integer :: log_unit = 6                   !< Unit for diagnostic messages (stdout)
    logical :: record_active = .false.        !< True if within a logical record
    logical :: buffer_valid = .false.         !< True if buffer contains unread data
    integer :: buffer_pos = 0                 !< Current read position in buffer
    integer :: buffer_len = 0                 !< Number of valid words in buffer
    integer(int32), allocatable :: buffer(:)  !< Word buffer for block data
    logical :: silent = .false.               !< Suppress error messages (for EOF drain)
  contains
    procedure :: init => output2_init
    procedure :: read_block_header => output2_read_block_header
    procedure :: read_words => output2_read_words
    procedure :: skip_record => output2_skip_record
  end type output2_reader

contains

  !===========================================================================
  ! C BINDING WRAPPER
  !===========================================================================

  !---------------------------------------------------------------------------
  !> Main entry point for C/C++ callers
  !>
  !> Converts Nastran modal data to MBDyn format. Input files (mbdyn.mat,
  !> mbdyn.tab) must exist in current directory.
  !>
  !> @param[in] outname                   Output filename (C string, null-terminated)
  !> @param[in] is_modal_displacement     Include initial modal displacements (0/1)
  !> @param[in] is_modal_velocity         Include initial modal velocities (0/1)
  !> @param[in] mass_direction            Lumped mass component selector (-1=auto, 0-3)
  !> @return                              0=success, non-zero=error code
  !---------------------------------------------------------------------------
  integer(c_int) function femgen_c_wrapper(outname,              &
                                          is_modal_displacement, &
                                          is_modal_velocity,     &
                                          mass_direction) bind(C, name='femgen')
    character(kind=c_char),    intent(in) :: outname(MAX_OUTNAME_LENGTH)
    integer(c_int32_t), value, intent(in) :: is_modal_displacement
    integer(c_int32_t), value, intent(in) :: is_modal_velocity
    integer(c_int32_t), value, intent(in) :: mass_direction
    
    type(femgen_model)                  :: model
    integer                             :: log_unit, idx, name_len, status
    integer                             :: opt_init_disp, opt_init_vel, opt_mass_idx
    character(len=MAX_OUTNAME_LENGTH-1) :: output_name
    character(len=:), allocatable       :: output_path

    femgen_c_wrapper = 0_c_int
    log_unit = 6
    output_name = ' '
    
    ! Convert C string to Fortran
    do idx = 1, min(len(output_name), size(outname))
      if (outname(idx) == c_null_char) exit
      output_name(idx:idx) = outname(idx)
    end do

    write(log_unit, '(A)') 'FEMGEN: translation of Binary NASTRAN modal data'
    write(log_unit, '(A)') '        for use in MBDyn''s modal joint element'
    write(log_unit, '(A)') ' '
    
    name_len = len_trim(output_name)
    if (name_len <= 0) then
      write(log_unit, '(A)') 'No output name provided'
      femgen_c_wrapper = 1_c_int
      return
    end if

    ! Read Nastran input files
    call femgen_read(FILE_MBDYN_MAT, FILE_MBDYN_TAB, model, log_unit, status)
    if (status /= 0) then
      femgen_c_wrapper = int(status, c_int)
      return
    end if

    ! Write MBDyn output file
    output_path = trim(output_name(1:name_len))
    if (.not. ends_with_fem(output_path)) output_path = output_path // '.fem'

    write(log_unit, '(A,A,A)') 'Output to file ''', output_path, ''''
    opt_init_disp = int(is_modal_displacement)
    opt_init_vel = int(is_modal_velocity)
    opt_mass_idx = int(mass_direction)
    call femgen_write(output_path, model, opt_init_disp, opt_init_vel, opt_mass_idx, status)
    if (status /= 0) then
      femgen_c_wrapper = int(status, c_int)
      return
    end if
  end function femgen_c_wrapper 

  !===========================================================================
  ! MAIN FILE READER
  !===========================================================================

  !---------------------------------------------------------------------------
  !> Read modal data from Nastran OUTPUT2/OUTPUT4 files
  !>
  !> Parses the binary Nastran files and populates the femgen_model structure
  !> with matrices (MHH, KHH, PHIX, LUMPMS) and tables (GPL, BGPDT).
  !>
  !> @param[in]    mat_path   Path to OUTPUT4 file (matrices)
  !> @param[in]    tab_path   Path to OUTPUT2 file (tables)
  !> @param[inout] model      Model structure to populate
  !> @param[in]    log_unit   Fortran unit for diagnostic output
  !> @param[out]   status     0=success, non-zero=error
  !---------------------------------------------------------------------------
  subroutine femgen_read(mat_path, tab_path, model, log_unit, status)
    character(len=*),   intent(in)    :: mat_path
    character(len=*),   intent(in)    :: tab_path
    type(femgen_model), intent(inout) :: model
    integer,            intent(in)    :: log_unit
    integer,            intent(out)   :: status

    type(output2_reader) :: reader
    real(real64), allocatable :: temp_matrix(:,:), lump_matrix(:,:)
    real :: timestamp(7)
    integer(int32), allocatable :: word_buffer(:)
    integer(int32) :: file_label(2), table_name(2)
    character(len=8) :: matrix_name
    
    integer :: out2_unit, out4_unit, io_ret, num_modes, words_read, ret_code
    integer :: node_idx, mode_idx, row_count, col_count, record_stride
    integer :: lump_rows, lump_cols, read_status
    logical :: out2_is_open, out4_is_open

    status = 0
    out2_is_open = .false.
    out4_is_open = .false.
    out2_unit = 11
    out4_unit = 12

    ! Open OUTPUT4 (matrix) file
    open(out4_unit, file=mat_path, status='old', form='unformatted', iostat=io_ret)
    if (io_ret /= 0) then
      write(log_unit,*) 'FEMGEN error while opening input file ', trim(mat_path), ', IRET=', io_ret
      call read_fail_and_cleanup(2, 'FEMGEN failed to open OUTPUT4 file')
      return
    end if
    out4_is_open = .true.

    ! Open OUTPUT2 (table) file
    open(out2_unit, file=tab_path, status='old', form='unformatted', iostat=io_ret)
    if (io_ret /= 0) then
      write(log_unit,*) 'FEMGEN error while opening input file ', trim(tab_path), ', IRET=', io_ret
      call read_fail_and_cleanup(2, 'FEMGEN failed to open OUTPUT2 file')
      return
    end if
    out2_is_open = .true.

    ! Read MHH (modal mass matrix)
    call output4_read_matrix(out4_unit, matrix_name, temp_matrix, row_count, col_count, read_status)
    write(log_unit,*) 'Reading Matrix ''', matrix_name, ''' ncol=', col_count, ' nrow=', row_count
    if (read_status /= 0) then
      call read_fail_and_cleanup(3, 'Error reading OUTPUT4 matrix')
      return
    end if
    if (matrix_name(1:3) /= 'MHH') then
      call read_fail_and_cleanup(3, 'Unexpected OUTPUT4 matrix order')
      return
    end if

    num_modes = row_count
    model%n_modes = num_modes
    allocate(model%mhh(num_modes, num_modes))
    model%mhh = temp_matrix(1:num_modes, 1:num_modes)
    deallocate(temp_matrix)

    ! Read KHH (modal stiffness matrix)
    call output4_read_matrix(out4_unit, matrix_name, temp_matrix, row_count, col_count, read_status)
    write(log_unit,*) 'Reading Matrix ''', matrix_name, ''' ncol=', col_count, ' nrow=', row_count
    if (read_status /= 0) then
      call read_fail_and_cleanup(3, 'Error reading OUTPUT4 matrix')
      return
    end if
    if (matrix_name(1:3) /= 'KHH') then
      call read_fail_and_cleanup(3, 'Unexpected OUTPUT4 matrix order')
      return
    end if

    allocate(model%khh(num_modes, num_modes))
    model%khh = 0.0_real64
    ! Handle diagonal eigenvalue format (ncol=5)
    if (row_count == num_modes .and. col_count == 5) then
      do node_idx = 1, num_modes
        do mode_idx = 1, num_modes
          model%khh(mode_idx, node_idx) = model%mhh(mode_idx, node_idx) * temp_matrix(node_idx, 1)
        end do
      end do
      ! Symmetrize
      do mode_idx = 1, num_modes
        do node_idx = mode_idx + 1, num_modes
          model%khh(mode_idx, node_idx) = 0.5_real64 * (model%khh(mode_idx, node_idx) + model%khh(node_idx, mode_idx))
          model%khh(node_idx, mode_idx) = model%khh(mode_idx, node_idx)
        end do
      end do
    else
      model%khh(1:row_count, 1:col_count) = temp_matrix(1:row_count, 1:col_count)
    end if
    deallocate(temp_matrix)

    ! Read LUMPMS (lumped mass diagonal)
    call output4_read_matrix(out4_unit, matrix_name, lump_matrix, lump_rows, lump_cols, read_status)
    write(log_unit,*) 'Reading Matrix ''', matrix_name, ''' ncol=', lump_cols, ' nrow=', lump_rows
    if (read_status /= 0) then
      call read_fail_and_cleanup(3, 'Error reading OUTPUT4 matrix')
      return
    end if
    if (matrix_name(1:6) /= 'LUMPMS') then
      call read_fail_and_cleanup(3, 'Unexpected OUTPUT4 matrix order')
      return
    end if

    allocate(model%lump_vec(lump_rows))
    model%lump_vec = lump_matrix(:, 1)
    deallocate(lump_matrix)

    ! Prepare word buffer for OUTPUT2 reading
    allocate(word_buffer(12*MAX_NODES))
    word_buffer = 0_int32

    ! Initialize OUTPUT2 reader and read first table
    call reader%init(out2_unit, log_unit, file_label, read_status)
    if (read_status /= 0) then
      call read_fail_and_cleanup(3, 'Error initializing OUTPUT2 reader')
      return
    end if
    call reader%read_block_header(table_name, timestamp, read_status)
    if (read_status /= 0) then
      call read_fail_and_cleanup(3, 'Error reading OUTPUT2 table header')
      return
    end if
    write(log_unit, '(A15,A4,A4)') ' Reading Table ', table_name(1), table_name(2)

    ! Check if we have GPL or OUG1 - require GPL/BGPDT for proper node data
    if (table_name(1) == transfer('OUG1', 0_int32) .or. &
        table_name(1) == transfer('OUGV', 0_int32)) then
      ! File starts with OUG1/OUGV - GPL/BGPDT missing
      call read_fail_and_cleanup(3, 'Error: GPL/BGPDT tables missing from OUTPUT2 file - check NASTRAN ALTER')
      return
    end if

    ! Read GPL table (Grid Point List)
      call reader%read_words(word_buffer, size(word_buffer), 0, words_read, ret_code)
      if (ret_code == -99) then
        call read_fail_and_cleanup(3, 'Error while reading GPL table')
        return
      else if (ret_code /= 1) then
        call read_fail_and_cleanup(3, 'Malformed GPL table trailer')
        return
      end if

      ! Extract node IDs from GPL table
      model%n_nodes = words_read
      allocate(model%node_ids(model%n_nodes))
      model%node_ids = word_buffer(1:model%n_nodes)

      allocate(model%xg(model%n_nodes), model%yg(model%n_nodes), model%zg(model%n_nodes))
      call reader%skip_record()
      call reader%skip_record()

      ! Read BGPDT table (grid point coordinates)
      call reader%read_block_header(table_name, timestamp, read_status)
      if (read_status /= 0) then
        call read_fail_and_cleanup(3, 'Error reading OUTPUT2 table header')
        return
      end if
      write(log_unit, '(A15,A4,A4)') ' Reading Table ', table_name(1), table_name(2)

      call reader%read_words(word_buffer, size(word_buffer), 0, words_read, ret_code)
      if (ret_code == -99) then
        call read_fail_and_cleanup(3, 'Error while reading BGPDT table')
        return
      else if (ret_code /= 1) then
        call read_fail_and_cleanup(3, 'Malformed BGPDT table trailer')
        return
      end if

      ! Determine BGPDT record format (12, 9, or 4 words per node)
      if (mod(words_read, 12) == 0) then
        record_stride = 12
      else if (mod(words_read, 9) == 0) then
        record_stride = 9
      else if (mod(words_read, 4) == 0) then
        record_stride = 4
      else
        write(log_unit,*) 'Unsupported BGPDT record layout', words_read
      call read_fail_and_cleanup(3, 'Unsupported BGPDT record layout')
      return
    end if

    ! Allocate internal IDs array
    allocate(model%internal_ids(model%n_nodes))

    ! Extract coordinates and internal IDs based on record format
    if (record_stride == 12) then
      do node_idx = 1, model%n_nodes
        call bgpdt12_get_xyz(word_buffer, node_idx, &
            model%xg(node_idx), model%yg(node_idx), model%zg(node_idx))
        ! BGPDT12: word 2 is internal ID
        model%internal_ids(node_idx) = word_buffer(12*(node_idx-1) + 2)
      end do
    else if (record_stride == 9) then
      do node_idx = 1, model%n_nodes
        call bgpdt9_get_xyz(word_buffer, node_idx, &
            model%xg(node_idx), model%yg(node_idx), model%zg(node_idx))
        ! BGPDT9: word 1 is internal ID (assuming similar structure)
        model%internal_ids(node_idx) = word_buffer(9*(node_idx-1) + 1)
      end do
    else
      do node_idx = 1, model%n_nodes
        call bgpdt4_get_xyz(word_buffer, node_idx, &
            model%xg(node_idx), model%yg(node_idx), model%zg(node_idx))
        ! BGPDT4: word 1 is internal ID
        model%internal_ids(node_idx) = word_buffer(4*(node_idx-1) + 1)
      end do
    end if

    ! Skip remaining BGPDT records
    call reader%skip_record()
    call reader%skip_record()

    ! Read OUG1 table (full eigenvectors for all nodes including dependent DOFs)
    ! OUG1 is REQUIRED - it provides complete mode shape data for all nodes
    call reader%read_block_header(table_name, timestamp, read_status)
    ! OUG1 in little-endian: 0x3147554F (ASCII 'OUG1' read as LE int32)
    if (read_status /= 0 .or. table_name(1) /= int(z'3147554F', int32)) then
      call read_fail_and_cleanup(3, 'Error: OUG1 table missing from OUTPUT2 file - check NASTRAN ALTER')
      return
    end if
    write(log_unit, '(A15,A4,A4)') ' Reading Table ', table_name(1), table_name(2)
    
    ! Allocate mode shape array from OUG1 (all nodes, all modes)
    allocate(model%phix(6*model%n_nodes, model%n_modes))
    model%phix = 0.0_real64
    model%phix_nrow = 6 * model%n_nodes
    model%phix_ncol = model%n_modes
    
    ! OUG1 table structure per mode:
    ! 1. IDENT record - header info (ACODE, TCODE, mode number, frequency, etc.)
    ! 2. DATA record(s) - nodal displacements (may be split across multiple records)
    ! 3. EOT marker (-1)
    
    ! Read each mode's data
    do mode_idx = 1, model%n_modes
      ! Skip the IDENT record (header for this mode)
      call reader%read_words(word_buffer, size(word_buffer), 0, words_read, ret_code)
      if (ret_code == -99) then
        call read_fail_and_cleanup(3, 'Error reading OUG1 IDENT record')
        return
      end if
      
      ! Read DATA record with actual displacement values
      call reader%read_words(word_buffer, size(word_buffer), 0, words_read, ret_code)
      if (ret_code == -99) then
        call read_fail_and_cleanup(3, 'Error reading OUG1 DATA record')
        return
      end if
      
      ! Parse the displacement data
      call parse_oug1_mode(word_buffer, words_read, model, mode_idx)
      
      ! Skip to end of this mode's data block (may be more records or EOT)
      do while (ret_code /= 1 .and. ret_code /= -99)
        call reader%read_words(word_buffer, size(word_buffer), 0, words_read, ret_code)
      end do
    end do
    
    write(log_unit, '(A,I0,A)') ' Successfully read OUG1 with ', model%n_modes, ' modes'
    
    ! Skip trailer records after OUG1
    call reader%skip_record()
    call reader%skip_record()

    ! Drain any remaining OUTPUT2 blocks (suppress error messages at EOF)
    reader%silent = .true.
    do
      call reader%read_words(word_buffer, size(word_buffer), 0, words_read, ret_code)
      if (ret_code == -99) then
        ! EOF or read error - this is normal after consuming all tables
        exit
      end if
      if (ret_code == 2) exit
    end do
    call cleanup_local()
    return

  contains

    subroutine read_fail_and_cleanup(error_code, error_msg)
      integer, intent(in) :: error_code
      character(len=*), intent(in) :: error_msg
      status = error_code
      if (len_trim(error_msg) > 0) write(log_unit, '(A)') trim(error_msg)
      call cleanup_local()
    end subroutine read_fail_and_cleanup

    subroutine cleanup_local()
      if (allocated(temp_matrix)) deallocate(temp_matrix)
      if (allocated(lump_matrix)) deallocate(lump_matrix)
      if (allocated(word_buffer)) deallocate(word_buffer)
      if (out2_is_open) then; close(out2_unit); out2_is_open = .false.; end if
      if (out4_is_open) then; close(out4_unit); out4_is_open = .false.; end if
    end subroutine cleanup_local
    
  end subroutine femgen_read 

  !===========================================================================
  ! MAIN FILE WRITER
  !===========================================================================

  !---------------------------------------------------------------------------
  !> Write MBDyn modal data file (.fem)
  !>
  !> Generates a formatted text file containing all modal data in the format
  !> expected by MBDyn's modal joint element.
  !>
  !> @param[in]    output_path       Output file path
  !> @param[inout] model             Model data (may be modified for mass fix)
  !> @param[in]    write_init_disp   Include initial modal displacements (0/1)
  !> @param[in]    write_init_vel    Include initial modal velocities (0/1)
  !> @param[in]    mass_component    Lumped mass treatment option
  !> @param[out]   status            0=success, non-zero=error
  !---------------------------------------------------------------------------
  subroutine femgen_write(output_path, model, write_init_disp, write_init_vel, mass_component, status)
    character(len=*),   intent(in)    :: output_path
    type(femgen_model), intent(inout) :: model
    integer,            intent(in)    :: write_init_disp, write_init_vel, mass_component
    integer,            intent(out)   :: status
    
    integer :: file_unit, io_ret, node_idx, mode_idx, dof_idx, row_end
    real(real64) :: mass_value
    
    status = 0
    file_unit = 10

    open(file_unit, file=output_path, iostat=io_ret)
    if (io_ret /= 0) then
      write(*,*) 'FEMGEN error opening output file, IRET=', io_ret
      status = 2
      return
    end if

    ! Record Group 1: Header
    write(file_unit, '(A24)') '** MBDyn MODAL DATA FILE'
    write(file_unit, '(A18)') '** NODE SET "ALL" '
    write(file_unit, *) ' '; write(file_unit, *) ' '
    write(file_unit, '(A25)') '** RECORD GROUP 1, HEADER'
    write(file_unit, '(A31, A39)') '**   REVISION,  NODE,  NORMAL, ', 'ATTACHMENT, CONSTRAINT, REJECTED MODES.'
    write(file_unit, '(6X, A4, 5I10)') 'REV0', model%n_nodes, model%n_modes, 0, 0, 0
    
    ! Record Group 2: Node list
    write(file_unit, '(A2)') '**'
    write(file_unit, '(A34,A9)') '** RECORD GROUP 2, FINITE ELEMENT ', 'NODE LIST'
    do node_idx = 1, model%n_nodes, 6
      row_end = min(node_idx + 5, model%n_nodes)
      write(file_unit, '(6I10)') (int(model%node_ids(mode_idx)), mode_idx=node_idx, row_end)
    end do

    ! Record Group 3: Initial modal displacements (optional)
    if (write_init_disp /= 0) then
      write(file_unit, '(A2)') '**'
      write(file_unit, '(A33, A13)') '** RECORD GROUP 3, INITIAL MODAL ', 'DISPLACEMENTS'
      write(file_unit, '(500(1X,1PE17.10))') (0.0*mode_idx, mode_idx=1, model%n_modes)
    end if
    
    ! Record Group 4: Initial modal velocities (optional)
    if (write_init_vel /= 0) then
      write(file_unit, '(A2)') '**'
      write(file_unit, '(A32,A10)') '** RECORD GROUP 4, INITIAL MODAL ', 'VELOCITIES'
      write(file_unit, '(500(1X,1PE17.10))') (0.0*mode_idx, mode_idx=1, model%n_modes)
    end if

    ! Record Groups 5-7: Nodal coordinates
    write(file_unit, '(A2)') '**'
    write(file_unit, '(A38)') '** RECORD GROUP 5, NODAL X COORDINATES'
    do node_idx = 1, model%n_nodes
      write(file_unit, '(E17.10)') snap_zero(real(model%xg(node_idx), kind=real64))
    end do
    write(file_unit, '(A2)') '**'
    write(file_unit, '(A38)') '** RECORD GROUP 6, NODAL Y COORDINATES'
    do node_idx = 1, model%n_nodes
      write(file_unit, '(E17.10)') snap_zero(real(model%yg(node_idx), kind=real64))
    end do
    write(file_unit, '(A2)') '**'
    write(file_unit, '(A38)') '** RECORD GROUP 7, NODAL Z COORDINATES'
    do node_idx = 1, model%n_nodes
      write(file_unit, '(E17.10)') snap_zero(real(model%zg(node_idx), kind=real64))
    end do

    ! Record Group 8: Mode shapes
    ! Since OUG1 is required, phix always contains full data for all nodes
    write(file_unit, '(A2)') '**'
    write(file_unit, '(A30)') '** RECORD GROUP 8, MODE SHAPES'

    ! Verify PHIX dimensions (should be 6*n_nodes rows, n_modes cols from OUG1)
    if (model%phix_nrow /= 6*model%n_nodes .or. model%phix_ncol /= model%n_modes) then
      write(*,*) 'PHIX dimensions mismatch: expected', 6*model%n_nodes, 'x', model%n_modes, &
                 ' got', model%phix_nrow, 'x', model%phix_ncol
      status = 3
      close(file_unit)
      return
    end if

    ! Write all mode shapes (PHIX stored as rows=DOFs, cols=modes)
    do mode_idx = 1, model%n_modes
      write(file_unit, '(A26,I2)') '**    NORMAL MODE SHAPE # ', mode_idx
      do node_idx = 1, model%n_nodes
        write(file_unit, '(6(1X,1PE17.10))') &
            (snap_zero(model%phix(6*(node_idx-1)+dof_idx, mode_idx)), dof_idx=1, 6)
      end do
    end do

    ! Record Group 9: Modal mass matrix
    write(file_unit, '(A2)') '**'
    write(file_unit, '(A36)') '** RECORD GROUP 9, MODAL MASS MATRIX'
    do mode_idx = 1, model%n_modes
      write(file_unit, '(500(1X,1PE17.10))') (snap_zero(model%mhh(mode_idx, node_idx)), node_idx=1, model%n_modes)
    end do
    
    ! Record Group 10: Modal stiffness matrix
    write(file_unit, '(A2)') '**'
    write(file_unit, '(A42)') '** RECORD GROUP 10, MODAL STIFFNESS MATRIX'
    do mode_idx = 1, model%n_modes
      write(file_unit, '(500(1X,1PE17.10))') (snap_zero(model%khh(mode_idx, node_idx)), node_idx=1, model%n_modes)
    end do

    ! Process lumped mass based on mass_component option
    if (mass_component == -1) then
      ! Auto-detect: check consistency across translational DOFs
      do node_idx = 0, size(model%lump_vec) - 6, 6
        if (abs(model%lump_vec(node_idx+1) - model%lump_vec(node_idx+2)) > FEMGEN_ZERO_TOL .or. &
            abs(model%lump_vec(node_idx+1) - model%lump_vec(node_idx+3)) > FEMGEN_ZERO_TOL) then
          mass_value = model%lump_vec(node_idx+1)
          if (abs(mass_value) <= FEMGEN_ZERO_TOL) mass_value = model%lump_vec(node_idx+2)
          if (abs(mass_value) <= FEMGEN_ZERO_TOL) mass_value = model%lump_vec(node_idx+3)
          do dof_idx = 1, 3
            if (abs(model%lump_vec(node_idx+dof_idx)) >= FEMGEN_ZERO_TOL .and. &
                abs(model%lump_vec(node_idx+dof_idx) - mass_value) > FEMGEN_ZERO_TOL) then
              write(*,*) 'Inconsistent lumped mass values for node #', node_idx/6 + 1
              status = 3
              close(file_unit)
              return
            end if
          end do
          do dof_idx = 1, 3
            model%lump_vec(node_idx+dof_idx) = snap_zero(mass_value)
          end do
        end if
      end do
    else if (mass_component /= 0) then
      ! Use specified component (1, 2, or 3) for all translational DOFs
      do node_idx = 0, size(model%lump_vec) - 6, 6
        mass_value = model%lump_vec(node_idx + mass_component)
        do dof_idx = 1, 3
          model%lump_vec(node_idx+dof_idx) = snap_zero(mass_value)
        end do
      end do
    end if

    ! Record Group 11: Lumped mass diagonal
    write(file_unit, '(A2)') '**'
    write(file_unit, '(A38,A12)') '** RECORD GROUP 11, DIAGONAL OF LUMPED', ' MASS MATRIX'
    do node_idx = 0, size(model%lump_vec) - 1, 6
      write(file_unit, '(500(1X,1PE17.10))') (snap_zero(model%lump_vec(node_idx+dof_idx)), dof_idx=1, 6)
    end do
    write(file_unit, '(A14)') '** END OF FILE'
    close(file_unit)
  end subroutine femgen_write

  !===========================================================================
  ! OUTPUT4 MATRIX READER
  !===========================================================================

  !---------------------------------------------------------------------------
  !> Read a matrix from Nastran OUTPUT4 binary file
  !>
  !> Parses a single matrix block from an OUTPUT4 file, handling both dense
  !> and sparse column formats, single and double precision data.
  !>
  !> @param[in]  file_unit    Fortran unit number (file must be open)
  !> @param[out] matrix_name  8-character matrix name from header
  !> @param[out] matrix       Allocated matrix data (nrow x ncol)
  !> @param[out] nrow         Number of rows
  !> @param[out] ncol         Number of columns
  !> @param[out] error_code   0=success, 1=bad type, 2=bad dims, 3=overflow, 4=I/O error
  !---------------------------------------------------------------------------
  subroutine output4_read_matrix(file_unit, matrix_name, matrix, nrow, ncol, error_code)
    integer, intent(in) :: file_unit
    character(len=8), intent(out) :: matrix_name
    real(real64), allocatable, intent(out) :: matrix(:,:)
    integer, intent(out) :: nrow, ncol
    integer, intent(out) :: error_code

    integer(int32) :: form_code, data_type       !< Matrix format and precision type
    integer(int32) :: col_idx, start_row         !< Current column and starting row
    integer(int32) :: num_words                  !< Words in current record
    integer        :: io_status, col_num, word_idx
    integer(int32), allocatable :: word_buffer(:)

    error_code = 0
    allocate(word_buffer(0))

    ! Read matrix header: ncol, nrow, form, type, name
    read(file_unit, iostat=io_status) ncol, nrow, form_code, data_type, matrix_name
    if (io_status /= 0) then
      call fail_with_code(4)
      return
    end if

    ! Validate dimensions
    if (ncol <= 0 .or. nrow <= 0) then
      call fail_with_code(2)
      return
    end if

    ! Check data type (1=single, 2=double)
    if (data_type > 2) then
      rewind(file_unit)
      call fail_with_code(1)
      return
    end if

    allocate(matrix(nrow, ncol))
    matrix = 0.0_real64

    ! Read column data
    do col_num = 1, ncol
      read(file_unit, iostat=io_status) col_idx, start_row, num_words
      if (io_status /= 0) then
        call fail_with_code(4)
        return
      end if
      backspace(file_unit)

      if (num_words < 0 .or. num_words > OUTPUT4_MAX_RECORD_WORDS) then
        call fail_with_code(3)
        return
      end if

      call ensure_buffer_size(int(num_words))
      read(file_unit, iostat=io_status) col_idx, start_row, num_words, (word_buffer(word_idx), word_idx=1, num_words)
      if (io_status /= 0) then
        call fail_with_code(4)
        return
      end if

      if (col_idx > ncol) exit

      if (start_row == 0) then
        call fill_sparse_column(nrow, data_type, col_idx, word_buffer(1:num_words), matrix)
      else
        call fill_dense_column(nrow, data_type, col_idx, start_row, word_buffer(1:num_words), matrix)
      end if
    end do

    ! Check for trailer record
    read(file_unit, iostat=io_status) col_idx, start_row, num_words
    if (io_status == 0) then
      if (col_idx <= ncol) backspace(file_unit)
    end if

    call cleanup()

  contains
    !> Set error code and cleanup
    subroutine fail_with_code(code)
      integer, intent(in) :: code
      error_code = code
      call cleanup()
    end subroutine fail_with_code

    !> Deallocate temporary buffer
    subroutine cleanup()
      if (allocated(word_buffer)) deallocate(word_buffer)
    end subroutine cleanup

    !> Ensure buffer has capacity for given number of words
    subroutine ensure_buffer_size(needed)
      integer, intent(in) :: needed
      integer(int32), allocatable :: temp(:)
      integer :: new_size
      if (needed <= 0) return
      if (allocated(word_buffer)) then
        if (size(word_buffer) >= needed) return
      end if
      if (.not. allocated(word_buffer)) then
        new_size = needed
      else if (size(word_buffer) <= 0) then
        new_size = needed
      else
        new_size = max(needed, size(word_buffer) * 2)
      end if
      allocate(temp(new_size))
      call move_alloc(temp, word_buffer)
    end subroutine ensure_buffer_size

    !---------------------------------------------------------------------------
    !> Fill matrix column with dense (contiguous) data
    !>
    !> @param[in]    nrows      Total matrix rows
    !> @param[in]    ntype      Data type (1=single, 2=double)
    !> @param[in]    icol       Column index
    !> @param[in]    irow0      Starting row index
    !> @param[in]    w          Raw word data
    !> @param[inout] a          Target matrix
    !---------------------------------------------------------------------------
    subroutine fill_dense_column(nrows, ntype, icol, irow0, w, a)
      integer, intent(in) :: nrows
      integer(int32), intent(in) :: ntype
      integer(int32), intent(in) :: icol, irow0
      integer(int32), intent(in) :: w(:)
      real(real64), intent(inout) :: a(:,:)
      integer :: val_idx, row, num_values
      if (ntype == 1) then
        do val_idx = 1, size(w)
          row = int(irow0) + val_idx - 1
          if (row >= 1 .and. row <= nrows) then
            a(row, int(icol)) = real(word_to_real32(w(val_idx)), real64)
          end if
        end do
      else
        num_values = size(w) / 2
        do val_idx = 1, num_values
          row = int(irow0) + val_idx - 1
          if (row >= 1 .and. row <= nrows) then
            a(row, int(icol)) = words_to_real64(w(2*val_idx-1), w(2*val_idx))
          end if
        end do
      end if
    end subroutine fill_dense_column

    !---------------------------------------------------------------------------
    !> Fill matrix column with sparse (string-encoded) data
    !>
    !> Sparse format uses "strings" where each starts with a header word:
    !> bits 31-16 = length L, bits 15-0 = starting row. Data follows header.
    !>
    !> @param[in]    nrows      Total matrix rows
    !> @param[in]    ntype      Data type (1=single, 2=double)
    !> @param[in]    icol       Column index
    !> @param[in]    w          Raw word data with string headers
    !> @param[inout] a          Target matrix
    !---------------------------------------------------------------------------
    subroutine fill_sparse_column(nrows, ntype, icol, w, a)
      integer, intent(in) :: nrows
      integer(int32), intent(in) :: ntype
      integer(int32), intent(in) :: icol
      integer(int32), intent(in) :: w(:)
      real(real64), intent(inout) :: a(:,:)
      integer :: string_pos, val_idx, row
      integer(int32) :: string_len, start_row, data_words
      string_pos = 1
      do
        if (string_pos > size(w)) exit
        ! Decode string header: length in upper 16 bits, start row in lower 16
        string_len = iand(ishft(w(string_pos), -16), int(z'FFFF', int32))
        start_row = iand(w(string_pos), int(z'FFFF', int32))
        data_words = string_len - 1
        if (ntype == 1) then
          do val_idx = 1, int(data_words)
            row = int(start_row) + val_idx - 1
            if (row >= 1 .and. row <= nrows) then
              a(row, int(icol)) = real(word_to_real32(w(string_pos + val_idx)), real64)
            end if
          end do
        else
          do val_idx = 1, int(data_words), 2
            row = int(start_row) + (val_idx-1)/2
            if (row >= 1 .and. row <= nrows) then
              a(row, int(icol)) = words_to_real64(w(string_pos + val_idx), w(string_pos + val_idx + 1))
            end if
          end do
        end if
        string_pos = string_pos + int(string_len)
        if (string_pos >= size(w)) exit
      end do
    end subroutine fill_sparse_column
  end subroutine output4_read_matrix

  ! ---------------------------------------------------------------------------
  ! OUTPUT2 Reader Logic
  ! ---------------------------------------------------------------------------

  !---------------------------------------------------------------------------
  !> Initialize OUTPUT2 reader and parse file header
  !>
  !> Validates the OUTPUT2 file format by reading the standard header sequence
  !> (date, timestamp, label) and prepares for table reading.
  !>
  !> @param[in]  file_unit   Fortran unit number (file must be open)
  !> @param[in]  log_unit    Optional unit for diagnostic output (default: 6)
  !> @param[out] label       2-word file label from header
  !> @param[out] status      0=success, non-zero=error
  !---------------------------------------------------------------------------
  subroutine output2_init(this, file_unit, log_unit, label, status)
    class(output2_reader), intent(inout) :: this
    integer, intent(in) :: file_unit
    integer, intent(in), optional :: log_unit
    integer(int32), intent(out) :: label(2)
    integer, intent(out) :: status
    
    integer(int32) :: key, month, day, year, timestamp(7)
    integer :: idx, io_status
    logical :: has_full_header
    
    status = 0
    this%file_unit = file_unit
    if (present(log_unit)) this%log_unit = log_unit
    if (.not. allocated(this%buffer)) allocate(this%buffer(200000))
    rewind(this%file_unit)

    ! 1. Header Key
    read(this%file_unit, iostat=io_status) key
    call check_read_error(io_status)
    if (status /= 0) return
    
    ! Check if file has full header (key=3) or simplified header (key=2)
    has_full_header = (key == 3_int32)
    
    if (.not. has_full_header) then
      ! Simplified format: rewind and skip full header parsing
      rewind(this%file_unit)
      label(1) = 0_int32
      label(2) = 0_int32
      this%record_active = .false.
      this%buffer_valid = .false.
      this%buffer_pos = 0
      this%buffer_len = 0
      return
    end if

    ! 2. Date
    read(this%file_unit, iostat=io_status) month, day, year
    call check_read_error(io_status)
    if (status /= 0) return

    ! 3. Info Key
    read(this%file_unit, iostat=io_status) key
    call check_read_error(io_status)
    if (status /= 0) return
    call check_key_value(key, 7_int32)
    if (status /= 0) return

    ! 4. Timestamp
    read(this%file_unit, iostat=io_status) (timestamp(idx), idx=1, 7)
    call check_read_error(io_status)
    if (status /= 0) return

    ! 5. Key=2
    read(this%file_unit, iostat=io_status) key
    call check_read_error(io_status)
    if (status /= 0) return
    call check_key_value(key, 2_int32)
    if (status /= 0) return

    ! 6. Label
    read(this%file_unit, iostat=io_status) label
    call check_read_error(io_status)
    if (status /= 0) return

    ! 7. Key=-1
    read(this%file_unit, iostat=io_status) key
    call check_read_error(io_status)
    if (status /= 0) return
    call check_key_value(key, -1_int32)
    if (status /= 0) return

    ! 8. Key=0
    read(this%file_unit, iostat=io_status) key
    call check_read_error(io_status)
    if (status /= 0) return
    call check_key_value(key, 0_int32)
    if (status /= 0) return

    this%record_active = .false.
    this%buffer_valid = .false.
    this%buffer_pos = 0
    this%buffer_len = 0

  contains
    subroutine check_read_error(io_stat)
      integer, intent(in) :: io_stat
      if (io_stat /= 0) then
        call fail_output2(this, 'IOPEN READ ERROR', 0_int32, status)
      end if
    end subroutine check_read_error

    subroutine check_key_value(val, expected)
      integer(int32), intent(in) :: val, expected
      if (val /= expected) then
        call fail_output2(this, 'IOPEN BAD KEY', val, status)
      end if
    end subroutine check_key_value

  end subroutine output2_init

  !---------------------------------------------------------------------------
  !> Read OUTPUT2 table block header
  !>
  !> Parses the header sequence that precedes each data table, extracting
  !> the table name and timestamp information.
  !>
  !> @param[out] table_name  2-word table identifier (e.g., 'GPL ', 'BGPDT')
  !> @param[out] timestamp   7-element timestamp array
  !> @param[out] status      0=success, non-zero=error
  !---------------------------------------------------------------------------
  subroutine output2_read_block_header(this, table_name, timestamp, status)
    class(output2_reader), intent(inout) :: this
    integer(int32),         intent(out)  :: table_name(2)
    real,                   intent(out)  :: timestamp(7)
    integer,                intent(out)  :: status

    integer(int32) :: key, header_info(2)
    integer :: io_status

    status = 0

    ! 1. Start Marker
    read(this%file_unit, iostat=io_status) key
    call check_read_error(io_status)
    if (status /= 0) return
    call check_key_value(key, 2_int32)
    if (status /= 0) return

    ! 2. Table Name
    read(this%file_unit, iostat=io_status) table_name
    call check_read_error(io_status)
    if (status /= 0) return

    ! 3. Separator
    read(this%file_unit, iostat=io_status) key
    call check_read_error(io_status)
    if (status /= 0) return
    call check_key_value(key, -1_int32)
    if (status /= 0) return

    ! 4. Info Marker
    read(this%file_unit, iostat=io_status) key
    call check_read_error(io_status)
    if (status /= 0) return
    call check_key_value(key, 7_int32)
    if (status /= 0) return

    ! 5. Timestamp
    read(this%file_unit, iostat=io_status) timestamp
    call check_read_error(io_status)
    if (status /= 0) return

    ! 6. Separator (some formats use -2, others use 3 - try to detect)
    read(this%file_unit, iostat=io_status) key
    call check_read_error(io_status)
    if (status /= 0) return
    if (key /= -2_int32 .and. key /= 3_int32) then
      call fail_output2(this, 'IHEADR BAD KEY after timestamp', key, status)
      return
    end if

    ! 7. Block Index - path depends on previous key
    if (key == 3_int32) then
      ! When key=3, next comes -2, then 1, then 0
      read(this%file_unit, iostat=io_status) key
      call check_read_error(io_status)
      if (status /= 0) return
      call check_key_value(key, -2_int32)
      if (status /= 0) return
    end if
    
    ! Now read block index (happens in both paths)
    read(this%file_unit, iostat=io_status) key
    call check_read_error(io_status)
    if (status /= 0) return
    call check_key_value(key, 1_int32)
    if (status /= 0) return

    ! 8. Separator
    read(this%file_unit, iostat=io_status) key
    call check_read_error(io_status)
    if (status /= 0) return
    call check_key_value(key, 0_int32)
    if (status /= 0) return

    ! 9. Length Check (this will be the next key regardless of path taken)
    read(this%file_unit, iostat=io_status) key
    call check_read_error(io_status)
    if (status /= 0) return
    if (key < 2_int32) then
      call fail_output2(this, 'IHEADR BAD KEY', key, status)
      return
    end if

    ! 10. Header Data
    read(this%file_unit, iostat=io_status) header_info
    call check_read_error(io_status)
    if (status /= 0) return

    ! 11. End Marker
    read(this%file_unit, iostat=io_status) key
    call check_read_error(io_status)
    if (status /= 0) return
    call check_key_value(key, -3_int32)
    if (status /= 0) return

    this%record_active = .false.
    this%buffer_valid = .false.
    this%buffer_pos = 0
    this%buffer_len = 0

  contains
    subroutine check_read_error(io_stat)
      integer, intent(in) :: io_stat
      if (io_stat /= 0) then
        call fail_output2(this, 'IHEADR READ ERROR', 0_int32, status)
      end if
    end subroutine check_read_error

    subroutine check_key_value(val, expected)
      integer(int32), intent(in) :: val, expected
      if (val /= expected) then
        call fail_output2(this, 'IHEADR BAD KEY', val, status)
      end if
    end subroutine check_key_value

  end subroutine output2_read_block_header

  !---------------------------------------------------------------------------
  !> Read words from OUTPUT2 data record
  !>
  !> Reads requested number of 32-bit words from the current data record,
  !> handling multi-block records transparently. Supports skipping and
  !> partial reads.\n  !>\n  !> @param[out]    dest_array     Destination array for read words
  !> @param[in]     words_requested Number of words to read (negative=skip)
  !> @param[in]     finish_record  If non-zero, skip to end of record after read
  !> @param[out]    words_read     Actual number of words transferred
  !> @param[out]    return_code    0=more data, 1=end of record, 2=end of table, -99=error
  !---------------------------------------------------------------------------
  subroutine output2_read_words(this, dest_array, words_requested, finish_record, words_read, return_code)
    class(output2_reader), intent(inout) :: this
    integer(int32), intent(inout) :: dest_array(:)
    integer, intent(in) :: words_requested
    integer, intent(in) :: finish_record
    integer, intent(out) :: words_read, return_code
    
    integer :: words_remaining, dest_pos, available, idx, io_status
    integer(int32) :: key

    words_remaining = words_requested
    dest_pos = 0
    return_code = 0
    words_read = 0

    do
      ! Start new record if needed
      if (.not. this%record_active) then
        read(this%file_unit, iostat=io_status) key
        call check_io(io_status, 'IREAD READ ERROR')
        if (return_code /= 0) return
        
        ! Check for end-of-table marker (key=0) - indicates all tables consumed
        if (key == 0_int32) then
          words_read = 0
          return_code = 2  ! End of table
          return
        end if
        
        call check_key(key, 1_int32)
        if (return_code /= 0) return

        read(this%file_unit, iostat=io_status) key
        call check_io(io_status, 'IREAD READ ERROR')
        if (return_code /= 0) return
        
        call check_key(key, 0_int32)
        if (return_code /= 0) return

        this%record_active = .true.
      end if

      ! Refill buffer if empty
      if (.not. this%buffer_valid) then
        read(this%file_unit, iostat=io_status) key
        call check_io_recover(io_status, 'IREAD READ ERROR')
        if (return_code /= 0) return

        if (key < 0_int32) then
          words_read = dest_pos
          return_code = 1
          this%record_active = .false.
          return
        end if
        if (key == 0_int32) then
          words_read = dest_pos
          return_code = 2
          this%record_active = .false.
          return
        end if
        if (key > int(OUTPUT2_MAX_BLOCK_WORDS, int32)) then
          write(this%log_unit, '(A,1X,I0)') 'IREAD BLOCK TOO LARGE:', key
          return_code = -99
          this%record_active = .false.
          return
        end if
        call ensure_buffer_capacity(this, int(key))
        read(this%file_unit, iostat=io_status) this%buffer(1:int(key))
        call check_io_recover(io_status, 'IREAD READ ERROR')
        if (return_code /= 0) return

        this%buffer_len = int(key)
        this%buffer_pos = 0
        this%buffer_valid = .true.
      end if
      
      ! Process request
      if (words_remaining < 0) then
        call skip_within_record(this, -words_remaining)
        words_remaining = 0
      else if (words_remaining == 0) then
        exit
      else
        available = this%buffer_len - this%buffer_pos
        if (available >= words_remaining) then
          do idx = 1, words_remaining
            dest_array(dest_pos + idx) = this%buffer(this%buffer_pos + idx)
          end do
          dest_pos = dest_pos + words_remaining
          this%buffer_pos = this%buffer_pos + words_remaining
          words_remaining = 0
        else
          do idx = 1, available
            dest_array(dest_pos + idx) = this%buffer(this%buffer_pos + idx)
          end do
          dest_pos = dest_pos + available
          words_remaining = words_remaining - available
          this%buffer_valid = .false.
          this%buffer_pos = 0
          this%buffer_len = 0
          cycle
        end if
      end if
      if (words_remaining == 0) exit
    end do
    
    if (finish_record /= 0) then
      call skip_to_end_of_record(this)
      words_read = dest_pos
      return_code = 1
      this%record_active = .false.
      return
    end if
    
    words_read = dest_pos
    if (this%buffer_valid) then
      return_code = 0
    else
      return_code = 1
      this%record_active = .false.
    end if

  contains
    subroutine check_io(io_stat, msg)
      integer, intent(in) :: io_stat
      character(len=*), intent(in) :: msg
      if (io_stat /= 0) then
        if (.not. this%silent) write(this%log_unit, '(A)') msg
        return_code = -99
      end if
    end subroutine check_io

    subroutine check_io_recover(io_stat, msg)
      integer, intent(in) :: io_stat
      character(len=*), intent(in) :: msg
      if (io_stat /= 0) then
        if (.not. this%silent) write(this%log_unit, '(A)') msg
        return_code = -99
        this%record_active = .false.
      end if
    end subroutine check_io_recover

    subroutine check_key(val, expected)
      integer(int32), intent(in) :: val, expected
      if (val /= expected) then
        if (.not. this%silent) write(this%log_unit, '(A,1X,I0)') 'IREAD BAD KEY:', val
        return_code = -99
      end if
    end subroutine check_key

    subroutine skip_within_record(this, num_skip)
      class(output2_reader), intent(inout) :: this
      integer, intent(in) :: num_skip
      integer :: remaining, io_stat
      integer(int32) :: blk_key
      remaining = num_skip
      do while (remaining > 0)
        if (.not. this%buffer_valid) then
          read(this%file_unit, iostat=io_stat) blk_key
          call check_io_recover(io_stat, 'IREAD READ ERROR')
          if (return_code /= 0) then
            this%buffer_valid = .false.
            return
          end if
          
          if (blk_key <= 0_int32) then
            this%record_active = .false.
            this%buffer_valid = .false.
            return
          end if
          if (blk_key > int(OUTPUT2_MAX_BLOCK_WORDS, int32)) then
            write(this%log_unit, '(A,1X,I0)') 'IREAD BLOCK TOO LARGE:', blk_key
            return_code = -99
            this%record_active = .false.
            this%buffer_valid = .false.
            return
          end if
          call ensure_buffer_capacity(this, int(blk_key))
          read(this%file_unit, iostat=io_stat) this%buffer(1:int(blk_key))
          call check_io_recover(io_stat, 'IREAD READ ERROR')
          if (return_code /= 0) then
            this%buffer_valid = .false.
            return
          end if
          this%buffer_len = int(blk_key)
          this%buffer_pos = 0
          this%buffer_valid = .true.
        end if
        if (this%buffer_len - this%buffer_pos >= remaining) then
          this%buffer_pos = this%buffer_pos + remaining
          remaining = 0
        else
          remaining = remaining - (this%buffer_len - this%buffer_pos)
          this%buffer_valid = .false.
          this%buffer_pos = 0
          this%buffer_len = 0
        end if
      end do
    end subroutine skip_within_record

    subroutine skip_to_end_of_record(this)
      class(output2_reader), intent(inout) :: this
      integer(int32) :: blk_key
      integer :: io_stat
      this%buffer_valid = .false.
      this%buffer_pos = 0
      this%buffer_len = 0
      do
        read(this%file_unit, iostat=io_stat) blk_key
        call check_io_recover(io_stat, 'IREAD READ ERROR')
        if (return_code /= 0) return

        if (blk_key <= 0_int32) exit
        if (blk_key > int(OUTPUT2_MAX_BLOCK_WORDS, int32)) then
          write(this%log_unit, '(A,1X,I0)') 'IREAD BLOCK TOO LARGE:', blk_key
          return_code = -99
          this%record_active = .false.
          return
        end if
        call ensure_buffer_capacity(this, int(blk_key))
        read(this%file_unit, iostat=io_stat) this%buffer(1:int(blk_key))
        call check_io_recover(io_stat, 'IREAD READ ERROR')
        if (return_code /= 0) return
      end do
    end subroutine skip_to_end_of_record
  end subroutine output2_read_words

  !---------------------------------------------------------------------------
  !> Skip current OUTPUT2 record
  !>\n  !> Advances file position past remaining data in current record.
  !---------------------------------------------------------------------------
  subroutine output2_skip_record(this)
    class(output2_reader), intent(inout) :: this
    integer(int32) :: dummy(1)
    integer :: num_read, ret_code
    ret_code = 0
    do while (ret_code == 0)
      call this%read_words(dummy, 0, 1, num_read, ret_code)
    end do
  end subroutine output2_skip_record

  !---------------------------------------------------------------------------
  !> Report OUTPUT2 read error
  !---------------------------------------------------------------------------
  subroutine fail_output2(this, msg, key, status)
    class(output2_reader), intent(in)  :: this
    character(len=*),      intent(in)  :: msg
    integer(int32),        intent(in)  :: key
    integer,               intent(out) :: status

    if (.not. this%silent) then
      if (key /= 0_int32) then
        write(this%log_unit, '(A,1X,I0)') trim(msg)//':', key
      else
        write(this%log_unit, '(A)') trim(msg)
      end if
    end if
    status = 1
  end subroutine fail_output2

  !---------------------------------------------------------------------------
  !> Ensure buffer has capacity for given number of words
  !---------------------------------------------------------------------------
  subroutine ensure_buffer_capacity(this, needed)
    class(output2_reader), intent(inout) :: this
    integer,               intent(in)    :: needed

    integer(int32), allocatable :: temp(:)
    integer :: current_size, new_size
    integer, parameter :: INITIAL_SIZE = 200000 

    if (needed <= 0) return
    if (.not. allocated(this%buffer)) then
      allocate(this%buffer(max(needed, INITIAL_SIZE)))
      return
    end if

    current_size = size(this%buffer)
    if (current_size >= needed) return
    new_size = max(needed, current_size * 2)
    allocate(temp(new_size))
    temp(1:current_size) = this%buffer
    call move_alloc(temp, this%buffer)
  end subroutine ensure_buffer_capacity

  !===========================================================================
  ! UTILITY FUNCTIONS
  !===========================================================================

  !---------------------------------------------------------------------------
  !> Convert 32-bit integer word to single-precision real
  !>
  !> Reinterprets the bit pattern of a 32-bit integer as IEEE single-precision
  !> floating point, used for reading OUTPUT4 single-precision data.
  !>
  !> @param[in] word  32-bit integer containing IEEE float bits
  !> @return          Equivalent single-precision real value
  !---------------------------------------------------------------------------
  pure function word_to_real32(word) result(x)
    integer(int32), intent(in) :: word
    real(real32) :: x
    x = transfer(word, x)
  end function word_to_real32

  !---------------------------------------------------------------------------
  !> Convert two 32-bit integer words to double-precision real
  !>
  !> Reinterprets a pair of 32-bit integers as IEEE double-precision
  !> floating point, used for reading OUTPUT4 double-precision data.
  !>
  !> @param[in] word1  Low 32 bits of IEEE double
  !> @param[in] word2  High 32 bits of IEEE double
  !> @return           Equivalent double-precision real value
  !---------------------------------------------------------------------------
  pure function words_to_real64(word1, word2) result(x)
    integer(int32), intent(in) :: word1, word2
    real(real64) :: x
    integer(int32) :: w(2)
    w(1) = word1
    w(2) = word2
    x = transfer(w, x)
  end function words_to_real64

  !---------------------------------------------------------------------------
  !> Snap very small values to exactly zero
  !>
  !> Values with magnitude below FEMGEN_ZERO_TOL are set to zero,
  !> improving output readability and avoiding numerical noise.
  !>
  !> @param[in] x  Input value
  !> @return       Zero if |x| <= tolerance, otherwise x unchanged
  !---------------------------------------------------------------------------
  pure elemental function snap_zero(x) result(y)
    real(real64), intent(in) :: x
    real(real64) :: y
    if (abs(x) <= FEMGEN_ZERO_TOL) then
      y = 0.0_real64
    else
      y = x
    end if
  end function snap_zero

  !---------------------------------------------------------------------------
  !> Parse OUG1 mode data and populate mode shape matrix
  !>
  !> OUG1 format: Each record contains eigenvector data for one mode.
  !> Data consists of sequences: [node_id, type_code, t1, t2, t3, r1, r2, r3]
  !> where node_id is int32, type_code is int32, and displacements are real32 or real64.
  !>
  !> @param[in]    words       Word buffer containing OUG1 record
  !> @param[in]    nwords      Number of words in buffer
  !> @param[inout] model       Model data structure to populate
  !> @param[in]    mode_idx    Current mode index (1-based)
  !---------------------------------------------------------------------------
  subroutine parse_oug1_mode(words, nwords, model, mode_idx)
    integer(int32), intent(in) :: words(:)
    integer, intent(in) :: nwords
    type(femgen_model), intent(inout) :: model
    integer, intent(in) :: mode_idx
    
    integer :: iword, node_id, node_idx, dof_base, words_per_node, found_idx
    real(real32) :: t1_32, t2_32, t3_32, r1_32, r2_32, r3_32
    real(real64) :: t1, t2, t3, r1, r2, r3
    
    ! Determine format: try 8 words/node (int32 + int32 + 6*real32)
    ! or 14 words/node (int32 + int32 + 6*real64)
    if (mod(nwords, 8) == 0) then
      words_per_node = 8  ! Real32 format
    else if (mod(nwords, 14) == 0) then
      words_per_node = 14 ! Real64 format
    else
      ! Try real32 anyway
      words_per_node = 8
    end if
    
    iword = 1
    do while (iword + words_per_node - 1 <= nwords)
      node_id = words(iword)  ! This is OUG1 point ID format: external_id * 10 + 3
      ! words(iword+1) is type code (typically 1 for GRID)
      
      ! Convert OUG1 point ID to external grid ID
      ! OUG1 uses format: point_id = external_id * 10 + device_code
      ! Device code is typically 1 or 3 for eigenvalue output
      ! Integer division node_id/10 gives external_id regardless of device code
      found_idx = 0
      do node_idx = 1, model%n_nodes
        ! Check if this OUG1 node_id matches this external node
        ! Formula: OUG1_id = external_id * 10 + device_code
        if (node_id / 10 == model%node_ids(node_idx)) then
          found_idx = node_idx
          exit
        end if
      end do
      
      if (found_idx > 0 .and. found_idx <= model%n_nodes) then
        dof_base = 6 * (found_idx - 1)
        
        if (words_per_node == 8) then
          ! Real32 format
          t1_32 = transfer(words(iword+2), t1_32)
          t2_32 = transfer(words(iword+3), t2_32)
          t3_32 = transfer(words(iword+4), t3_32)
          r1_32 = transfer(words(iword+5), r1_32)
          r2_32 = transfer(words(iword+6), r2_32)
          r3_32 = transfer(words(iword+7), r3_32)
          
          model%phix(dof_base+1, mode_idx) = real(t1_32, real64)
          model%phix(dof_base+2, mode_idx) = real(t2_32, real64)
          model%phix(dof_base+3, mode_idx) = real(t3_32, real64)
          model%phix(dof_base+4, mode_idx) = real(r1_32, real64)
          model%phix(dof_base+5, mode_idx) = real(r2_32, real64)
          model%phix(dof_base+6, mode_idx) = real(r3_32, real64)
        else
          ! Real64 format
          t1 = words_to_real64(words(iword+2), words(iword+3))
          t2 = words_to_real64(words(iword+4), words(iword+5))
          t3 = words_to_real64(words(iword+6), words(iword+7))
          r1 = words_to_real64(words(iword+8), words(iword+9))
          r2 = words_to_real64(words(iword+10), words(iword+11))
          r3 = words_to_real64(words(iword+12), words(iword+13))
          
          model%phix(dof_base+1, mode_idx) = t1
          model%phix(dof_base+2, mode_idx) = t2
          model%phix(dof_base+3, mode_idx) = t3
          model%phix(dof_base+4, mode_idx) = r1
          model%phix(dof_base+5, mode_idx) = r2
          model%phix(dof_base+6, mode_idx) = r3
        end if
      end if
      
      iword = iword + words_per_node
    end do
  end subroutine parse_oug1_mode

  !---------------------------------------------------------------------------
  !> Extract node IDs and coordinates from OUG1 eigenvector table
  !>
  !> When GPL/BGPDT tables are missing, this extracts node information
  !> directly from the first mode in the OUG1 displacement table.
  !>
  !> @param[inout] reader      OUTPUT2 reader object  
  !> @param[inout] word_buffer Work buffer for reading
  !> @param[inout] model       Modal model data structure to populate
  !> @param[in]    log_unit    Unit for diagnostic output
  !> @param[out]   status      0=success, non-zero=error
  !---------------------------------------------------------------------------
  subroutine extract_nodes_from_oug1(reader, word_buffer, model, log_unit, status)
    type(output2_reader), intent(inout) :: reader
    integer(int32), intent(inout) :: word_buffer(:)
    type(femgen_model), intent(inout) :: model
    integer, intent(in) :: log_unit
    integer, intent(out) :: status
    
    integer :: words_read, ret_code, node_idx, iword, words_per_node
    integer(int32) :: node_id, type_code
    integer :: max_nodes, n_nodes_found
    integer(int32), allocatable :: temp_node_ids(:)
    real(real64), allocatable :: temp_xg(:), temp_yg(:), temp_zg(:)
    
    status = 0
    max_nodes = size(word_buffer) / 8  ! Conservative estimate
    allocate(temp_node_ids(max_nodes))
    allocate(temp_xg(max_nodes), temp_yg(max_nodes), temp_zg(max_nodes))
    
    ! Read first mode from OUG1 to extract node IDs
    ! OUG1 format: each node has [node_id, type_code, t1, t2, t3, r1, r2, r3]
    ! All values are 32-bit (integer for IDs, real32 for displacements)
    call reader%read_words(word_buffer, size(word_buffer), 0, words_read, ret_code)
    if (ret_code == -99) then
      status = 1
      write(log_unit, '(A)') ' Error reading OUG1 data for node extraction'
      return
    end if
    
    ! Parse OUG1 records - typically 8 words per node
    words_per_node = 8
    n_nodes_found = words_read / words_per_node
    
    if (n_nodes_found == 0) then
      status = 1
      write(log_unit, '(A)') ' No nodes found in OUG1 table'
      return
    end if
    
    ! Extract node IDs (coordinates will be zero - PHIX from OUTPUT4 doesn't need them)
    iword = 1
    do node_idx = 1, n_nodes_found
      if (iword + words_per_node - 1 > words_read) exit
      
      node_id = word_buffer(iword)
      type_code = word_buffer(iword + 1)
      
      ! Store node ID (coordinates set to zero since we use PHIX from OUTPUT4)
      temp_node_ids(node_idx) = node_id
      temp_xg(node_idx) = 0.0_real64
      temp_yg(node_idx) = 0.0_real64
      temp_zg(node_idx) = 0.0_real64
      
      iword = iword + words_per_node
    end do
    
    ! Allocate and copy to model
    model%n_nodes = n_nodes_found
    allocate(model%node_ids(model%n_nodes))
    allocate(model%xg(model%n_nodes), model%yg(model%n_nodes), model%zg(model%n_nodes))
    
    model%node_ids = temp_node_ids(1:n_nodes_found)
    model%xg = temp_xg(1:n_nodes_found)
    model%yg = temp_yg(1:n_nodes_found)
    model%zg = temp_zg(1:n_nodes_found)
    
    write(log_unit, '(A,I0,A)') ' Extracted ', n_nodes_found, ' nodes from OUG1'
    
    ! Skip remaining records for this mode
    if (ret_code == 1) then
      call reader%skip_record()
    end if
    
    ! Skip to end of OUG1 table
    do
      call reader%read_words(word_buffer, size(word_buffer), 0, words_read, ret_code)
      if (ret_code == -99) then
        exit  ! Error or end of file
      else if (ret_code == 1) then
        call reader%skip_record()  ! Skip trailer
      end if
      if (words_read == 0) exit
    end do
    
    deallocate(temp_node_ids, temp_xg, temp_yg, temp_zg)
  end subroutine extract_nodes_from_oug1

  !---------------------------------------------------------------------------
  !> Convert single character to lowercase
  !>
  !> @param[in] c  Single character (any case)
  !> @return       Lowercase equivalent (A-Z -> a-z), others unchanged
  !---------------------------------------------------------------------------
  pure function lower1(c) result(lc)
    character(len=1), intent(in) :: c
    character(len=1) :: lc
    integer :: ia
    ia = iachar(c)
    if (ia >= iachar('A') .and. ia <= iachar('Z')) then
      lc = achar(ia + (iachar('a') - iachar('A')))
    else
      lc = c
    end if
  end function lower1

  !---------------------------------------------------------------------------
  !> Check if filename ends with .fem extension (case-insensitive)
  !>
  !> @param[in] name  Filename string
  !> @return          True if name ends with .fem/.FEM/.Fem etc.
  !---------------------------------------------------------------------------
  pure function ends_with_fem(name) result(ok)
    character(len=*), intent(in) :: name
    logical :: ok
    integer :: n
    ok = .false.
    n = len_trim(name)
    if (n < 4) return
    ok = (lower1(name(n-3:n-3)) == '.') .and. &
         (lower1(name(n-2:n-2)) == 'f') .and. &
         (lower1(name(n-1:n-1)) == 'e') .and. &
         (lower1(name(n:n))     == 'm')
  end function ends_with_fem 

  !---------------------------------------------------------------------------
  !> Extract XYZ coordinates from 12-word BGPDT record (double precision)
  !>
  !> BGPDT format with 12 words per node: ID(1), CP(1), CD(1), PS(1), ICS(1),
  !> reserved(1), X(2), Y(2), Z(2). Coordinates are 64-bit (2 words each).
  !>
  !> @param[in]  words  Word array containing BGPDT data
  !> @param[in]  idx    1-based node index
  !> @param[out] x,y,z  Global coordinates
  !---------------------------------------------------------------------------
  subroutine bgpdt12_get_xyz(words, idx, x, y, z)
    integer(int32), intent(in) :: words(:)
    integer, intent(in) :: idx
    real(real64), intent(out) :: x, y, z
    integer :: base
    base = 12*(idx-1) + 1
    x = words_to_real64(words(base+6), words(base+7))
    y = words_to_real64(words(base+8), words(base+9))
    z = words_to_real64(words(base+10), words(base+11))
  end subroutine bgpdt12_get_xyz

  !---------------------------------------------------------------------------
  !> Extract XYZ coordinates from 9-word BGPDT record (single precision)
  !>
  !> BGPDT format with 9 words per node: ID(1), CP(1), CD(1), PS(1), ICS(1),
  !> reserved(1), X(1), Y(1), Z(1). Coordinates are 32-bit (1 word each).
  !>
  !> @param[in]  words  Word array containing BGPDT data
  !> @param[in]  idx    1-based node index
  !> @param[out] x,y,z  Global coordinates
  !---------------------------------------------------------------------------
  subroutine bgpdt9_get_xyz(words, idx, x, y, z)
    integer(int32), intent(in) :: words(:)
    integer, intent(in) :: idx
    real(real64), intent(out) :: x, y, z
    integer :: base
    base = 9*(idx-1) + 1
    x = real(word_to_real32(words(base+6)), real64)
    y = real(word_to_real32(words(base+7)), real64)
    z = real(word_to_real32(words(base+8)), real64)
  end subroutine bgpdt9_get_xyz

  !---------------------------------------------------------------------------
  !> Extract XYZ coordinates from 4-word BGPDT68 record (single precision)
  !>
  !> Compact BGPDT68 format with 4 words per node: ID(1), X(1), Y(1), Z(1).
  !> Coordinates are 32-bit (1 word each).
  !>
  !> @param[in]  words  Word array containing BGPDT68 data
  !> @param[in]  idx    1-based node index
  !> @param[out] x,y,z  Global coordinates
  !---------------------------------------------------------------------------
  subroutine bgpdt4_get_xyz(words, idx, x, y, z)
    integer(int32), intent(in) :: words(:)
    integer, intent(in) :: idx
    real(real64), intent(out) :: x, y, z
    integer :: base
    base = 4*(idx-1) + 1
    x = real(word_to_real32(words(base+1)), real64)
    y = real(word_to_real32(words(base+2)), real64)
    z = real(word_to_real32(words(base+3)), real64)
  end subroutine bgpdt4_get_xyz

end module femgen_mod
