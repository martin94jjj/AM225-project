# Load the common configuration file
include config.mk
#objs=poisson_fft.o file_output.o
src=$(patsubst %.o,%.cc,$(objs))
execs=cv cv_dx_conv cv_example cv_example_bck cv_example_cn pfft_test 

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

#pfft_test: pfft_test.cc $(objs)
	#$(cxx) $(cflags) $(fftw_iflags) -o $@ $^ $(fftw_lflags) $(lp_lflags)

clean:
	rm -f $(execs) $(objs) libf2d.a

.PHONY: clean depend
