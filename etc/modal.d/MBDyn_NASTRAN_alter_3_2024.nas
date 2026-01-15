$ $Header$
$ Modified for MSC/NASTRAN 2024.1  
$ SOL 103 MODAL ANALYSIS - Output files for MBDyn
$
$ MBDyn (C) is a multibody analysis code.
$ http://www.mbdyn.org
$ 
$ Copyright (C) 1996-2023 Pierangelo Masarati
$
$ Pierangelo Masarati     <pierangelo.masarati@polimi.it>
$ Paolo Mantegazza        <paolo.mantegazza@polimi.it>
$
$ Author: Giuseppe Quaranta <giuseppe.quaranta@polimi.it>
$         Alessandro Cocco <alessandro.cocco@polimi.it>
$         Dipartimento di Ingegneria Aerospaziale - Politecnico di Milano
$
$ Updated for MSC Nastran 2024.1
$ Using correct datablock names from 2024.1 DMAP structure
$
$  ASSIGN FILES FOR OUTPUT2 AND OUTPUT4
$
ASSIGN OUTPUT2='mbdyn.tab' STATUS=UNKNOWN UNIT=11  $ tables
ASSIGN OUTPUT4='mbdyn.mat' STATUS=UNKNOWN UNIT=15  $ matrices  
$
$$$ EXECUTIVE CONTROL DECK
$
TIME 500
SOL 103
$
$--------------------------------------------------------------------------
$ Strategy based on 2024.1 DMAP structure:
$ - GPL, BGPDT from PHASE0
$ - MIX (modal mass) and EIGVMAT (eigenvalues) from XREAD after READ module  
$ - OUGV1 (output eigenvectors) from SEDRCVR
$
$ MIX = Modal mass matrix (nmodes x nmodes)
$ KHH = Modal stiffness/generalized stiffness (nmodes x nmodes), derived from LAMA via LAMX
$ LUMPMS will be diagonal of MGG
$--------------------------------------------------------------------------
$
$ Output geometry tables
$ Write GPL and BGPDT during PHASE0 (geometry stage)
COMPILE      PHASE0
ALTER        615
OUTPUT2      GPL,BGPDT,,,//-1/11  $
$
$ Output geometry and modal matrices after READ module computes them
COMPILE      XREAD
ALTER        131
$
$ Rename MIX to MHH (modal mass) for femgen compatibility
EQUIVX       MIX/MHH/-1 $

$ Build KHH (diagonal eigenvalues matrix) from LAMA table
LAMX         , ,LAMA/KHH/-1 $
$
$ Output MHH and KHH (leave OUTPUT4 open; LUMPMS appended later)
OUTPUT4      MHH,KHH,,//-1/15 $
ENDALTER

$ Compute and append LUMPMS in a later subdmap where the physical mass matrix exists
COMPILE      SEDRCVR
ALTER        2700
DIAGONAL     MGG/LUMPMS/'COLUMN'/1.  $
OUTPUT4      LUMPMS,,,,//-2/15 $
OUTPUT4      ,,,//-9/15 $

$ Output full eigenvectors (OUG1) with dependent DOFs recovered
OUTPUT2      OUG1,,,//-2/11 $
OUTPUT2      ,,,,//-9/11 $
ENDALTER
