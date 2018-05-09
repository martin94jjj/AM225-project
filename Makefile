# Load the common configuration file
include config.mk
#objs=poisson_fft.o file_output.o
src=$(patsubst %.o,%.cc,$(objs))
execs=cv cv_dx_conv cv_example cv_example_bck cv_example_cn 
execs=cv_example_bck_conv cv_example_cn_conv cv_example_bck_changing_x 
execs=cv_example_bck_changing_x_conv
execs=butler_volmer butler_volmer_changing_alpha cv_example_bck_changing_x_peak_separation

execs=butler_volmer_simulate_NQ butler_volmer_find_peak_potential

all: $(execs)


cv: cv.cc 
	$(cxx) $(cflags) $(iflags) -o $@ $^ $(lp_lflags)

cv_dx_conv: cv_dx_conv.cc
	$(cxx) $(cflags) $(iflags) -o $@ $^ $(lp_lflags)

cv_example: cv_example.cc 
	$(cxx) $(cflags) $(iflags) -o $@ $^ $(lp_lflags)

cv_example_bck: cv_example_bck.cc
	$(cxx) $(cflags) $(iflags) -o $@ $^ $(lp_lflags)

cv_example_cn: cv_example_cn.cc
	$(cxx) $(cflags) $(iflags) -o $@ $^ $(lp_lflags)

cv_example_bck_conv: cv_example_bck_conv.cc
	$(cxx) $(cflags) $(iflags) -o $@ $^ $(lp_lflags)

cv_example_cn_conv: cv_example_cn_conv.cc
	$(cxx) $(cflags) $(iflags) -o $@ $^ $(lp_lflags)

cv_example_bck_changing_x: cv_example_bck_changing_x.cc
	$(cxx) $(cflags) $(iflags) -o $@ $^ $(lp_lflags)

cv_example_bck_changing_x_conv: cv_example_bck_changing_x_conv.cc
	$(cxx) $(cflags) $(iflags) -o $@ $^ $(lp_lflags)

butler_volmer: butler_volmer.cc
	$(cxx) $(cflags) $(iflags) -o $@ $^ $(lp_lflags)

butler_volmer_changing_alpha: butler_volmer_changing_alpha.cc
	$(cxx) $(cflags) $(iflags) -o $@ $^ $(lp_lflags)

cv_example_bck_changing_x_peak_separation: cv_example_bck_changing_x_peak_separation.cc
	$(cxx) $(cflags) $(iflags) -o $@ $^ $(lp_lflags)

butler_volmer_simulate_NQ: butler_volmer_simulate_NQ.cc
	$(cxx) $(cflags) $(iflags) -o $@ $^ $(lp_lflags)

butler_volmer_find_peak_potential: butler_volmer_find_peak_potential.cc
	$(cxx) $(cflags) $(iflags) -o $@ $^ $(lp_lflags)

clean:
	rm -f $(execs) $(objs) libf2d.a

.PHONY: clean depend
