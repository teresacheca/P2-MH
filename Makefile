SRCPATH := src
OBJPATH := obj

result := pracs_exe
OBJS := $(shell cd $(SRCPATH) && find . -type f -iname '*.cpp' | sed 's/\.cpp/\.o/g' | xargs -I {} echo "$(OBJPATH)/"{})

uname := $(shell uname -s)
en_macos := $(findstring Darwin,$(uname))
en_linux := $(findstring Linux,$(uname))
CXX := $(if $(en_linux), g++, clang++ )
sistoper := $(if $(en_macos), macOS, Linux )

cc_flags_common := -std=c++11 -O3 -I/usr/include -I$(SRCPATH)/
cc_flags_linux := -DLINUX
cc_flags_macos := -DMACOS
CXXFLAGS := $(cc_flags_common) $(if $(en_linux), $(cc_flags_linux), $(cc_flags_macos))

glu_flag_macos := /System/Library/Frameworks/OpenGL.framework/Versions/A/Libraries/libGLU.dylib
ld_libs_common := -ljpeg
ld_libs_linux := -lGLEW -lGLU -lglut -lGL
ld_libs_macos := -framework OpenGL -framework GLUT $(glu_flag_macos)
ld_libs := $(ld_libs_common) $(if $(en_linux), $(ld_libs_linux), $(ld_libs_macos))

.PHONY: run dev clean mrproper
run: $(result)
	./$(result)

dev: $(result)
	./$(result)

$(result): $(OBJPATH) $(OBJS)
	$(CXX) -o $(result) $(OBJS) $(ld_libs)

-include $(OBJS:.o=.o.d)

$(OBJPATH)/%.o: $(SRCPATH)/%.cpp
	$(CXX) -c -o $@ $< $(CXXFLAGS)
	$(CXX) -MM $< -o $@.d.tmp $(CXXFLAGS)
	@sed -e 's|.*:|$@:|' < $@.d.tmp > $@.d
	@rm $@.d.tmp

$(OBJPATH):
	cd $(SRCPATH) && find . -type d -exec mkdir -p ../$(OBJPATH)/{} \;

clean:
	rm -rf $(OBJPATH) $(result)

mrproper: clean
