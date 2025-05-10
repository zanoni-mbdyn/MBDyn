try
  clear all;
  close all;

  pkg_prefix = getenv("OCT_PKG_INSTALL_PREFIX");

  if (~isempty(pkg_prefix))
    pkg("local_list", fullfile(pkg_prefix, "octave_packages"));
  endif

  pkg load mboct-mbdyn-pkg;

  ad = [false, true];
  static = [false, true];
  elastic = [0, 1, 2];

  for k=1:numel(ad)
    fd = -1;
    unwind_protect
      fd = fopen("beam_ad.cod", "wt");
      if (fd == -1)
        error("failed to open file");
      endif
      if (ad(k))
        fputs(fd, "use automatic differentiation;\n");
      else
        fputs(fd, "## use automatic differentiation;\n");
      endif
    unwind_protect_cleanup
      if (fd ~= -1)
        fclose(fd);
      endif
    end_unwind_protect
    for l=1:numel(static)
      fd = -1;
      unwind_protect
        fd = fopen("beam_static.cod", "wt");
        if (fd == -1)
          error("failed to open file");
        endif
        if (static(l))
          fputs(fd, "model: static;\n");
        else
          fputs(fd, "## model: static;\n");
        endif
      unwind_protect_cleanup
        fclose(fd);
      end_unwind_protect
      for m=1:numel(elastic)
        fd = -1;
        unwind_protect
          fd = fopen("beam_const_law.csl", "wt");
          if (fd == -1)
            error("failed to open file");
          endif
          switch (elastic(m))
            case 0
              fputs(fd, "include: \"beam_elastic.csl\";\n");
            case 1
              fputs(fd, "include: \"beam_viscoelastic.csl\";\n");
            case 2
              fputs(fd, "include: \"beam.csl\";\n");
          endswitch
        unwind_protect_cleanup
          fclose(fd);
        end_unwind_protect

        opts.verbose = false;
        opts.output_file = "stage1";
        opts.logfile = "stage1.stdout";
        mbdyn_solver_run("beam_stage1.mbd", opts);
        opts.output_file = "stage2";
        opts.logfile = "stage2.stdout";
        mbdyn_solver_run("beam_stage2.mbd", opts);

        [t1, traj1, def1, vel1, acc1, node_id1] = mbdyn_post_load_output_struct("stage1");
        [t2, traj2, def2, vel2, acc2, node_id2] = mbdyn_post_load_output_struct("stage2");

        max_diff = 0;

        for j=1:numel(traj1)
          x2 = traj2{j};
          x1 = zeros(size(x2));

          for i=1:columns(x1)
            x1(:, i) = interp1(t1, traj1{j}(:, i), t2);
          endfor

          curr_diff = max(max(abs(x1 - x2))) / max(1, max(max(abs(x1))));
          max_diff = max(max_diff, curr_diff);
        endfor

        fprintf(stderr, "[%d,%d,%d] maximum difference between stage1 and stage2: %e\n", k, l, m, max_diff);
        assert(max_diff < 1e-10);
      endfor
    endfor
  endfor
catch
  gtest_error = lasterror();
  gtest_fail(gtest_error, __FILE__);
  rethrow(gtest_error);
end_try_catch
