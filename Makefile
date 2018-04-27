# Load the common configuration file
include config.mk

src=$(patsubst %.o,%.cc,$(objs))
execs=cv cv_example

all: $(execs)


cv: cv.cc 
	$(cxx) $(cflags) $(iflags) -o $@ $^ $(lp_lflags)

cv_example: cv_example.cc 
	$(cxx) $(cflags) $(iflags) -o $@ $^ $(lp_lflags)


clean:
	rm -f $(execs) $(objs) libf2d.a

.PHONY: clean depend
