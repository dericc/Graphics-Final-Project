#
# Makefile for assignment #4
#

#
# Compile and link options
#

CXX=g++
CXXFLAGS=-Wall -I. -O3 -DUSE_JPEG



#
# OpenGL libraries
#
UNAME := $(shell uname)
ifneq (,$(findstring Darwin,$(UNAME)))
	GLLIBS = -framework GLUT -framework OpenGL
        CXXFLAGS := $(CXXFLAGS) -mmacosx-version-min=10.8
else
  ifneq (,$(findstring CYGWIN,$(UNAME)))
	GLLIBS = -lopengl32 -lglu32 -lwinmm -lgdi32
  else
	GLLIBS = -lGLU -lGL
  endif
endif



#
# GNU Make: targets that don't build files
#

.PHONY: all clean distclean

#
# Rules encoding targets and dependencies.  By default, the first of
# these is built, but you can also build any individual target by
# passing it to make - e.g., "make imgpro" or "make clean"
#
# Notice that many of the dependencies are implicit (e.g. a .o depends
# on its corresponding .cpp), as are many of the compilation rules.
#

OBJS=game.o particle.o R3Scene.o raytrace.o
LIBS=R3/libR3.a R2/libR2.a jpeg/libjpeg.a fglut/libfglut.a
ifneq (,$(findstring Linux,$(UNAME)))
	OPTS =  -I"../irrKlang/include" -L"../irrKlang/bin/linux-gcc-64" ../irrKlang/bin/linux-gcc-64/libIrrKlang.so
	INSTALL_NAME_TOOL = 
else 
	OPTS =  -I"../irrKlang/include" -L"../irrKlang/bin/macosx-gcc" ../irrKlang/bin/macosx-gcc/libirrklang.dylib
	INSTALL_NAME_TOOL = install_name_tool -change /usr/local/lib/libirrklang.dylib $(shell pwd)/../irrKlang/bin/macosx-gcc/libirrklang.dylib game
endif

all: game

R3/libR3.a: 
	    $(MAKE) -C R3

R2/libR2.a: 
	    $(MAKE) -C R2

jpeg/libjpeg.a: 
	    $(MAKE) -C jpeg

fglut/libfglut.a: 
	$(MAKE) -C fglut

game: $(OBJS) $(LIBS) 
	    rm -f $@ 
	    $(CXX) $(CXXFLAGS) $^ -lm -o $@ $(GLLIBS) $(OPTS)
	    $(INSTALL_NAME_TOOL)
	    
clean:
	    rm -f *.o game
		$(MAKE) -C R3 clean
		$(MAKE) -C R2 clean
		$(MAKE) -C jpeg clean
		$(MAKE) -C fglut clean


distclean:  clean

