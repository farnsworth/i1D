

CODE_DIR=code/

OBJECTS=${CODE_DIR}system.f90 \
	${CODE_DIR}wl.f90 \
	${CODE_DIR}exact.f90 \
	${CODE_DIR}quench.f90 \
	${CODE_DIR}canonical.f90 \
	${CODE_DIR}wl_new.f90

MODULE=ising1D

default :
	@echo 'to produce the module type at the shell prompt:'
	@echo 'make target'
	@echo 'Where target is one of the following:'
	@echo '  module       to produce the interface'
	@echo '  compile      to compile and produce the wrap'
	@echo '  all          to do the two things in a single step'
	@echo '  clean        to clean .pyf and .so files'
	@echo '  veryclean    same of clean plus the remove'
	@echo '               of pyc and temp files'

all:
	make module
	make compile

module:$(OBJECTS)
	rm -f $(CODE_DIR)$(MODULE).pyf
	f2py -m $(MODULE) -h $(CODE_DIR)$(MODULE).pyf $(OBJECTS)

compile:$(OBJECTS) $(CODE_DIR)$(MODULE).pyf
	f2py -c --fcompiler=gnu95 $(CODE_DIR)$(MODULE).pyf $(OBJECTS)

clean:
	rm $(CODE_DIR)$(MODULE).pyf
	rm $(MODULE).so

veryclean:
	rm -f $(CODE_DIR)$(MODULE).pyf
	rm -f $(MODULE).so
	rm -f *~
	rm -f */*~
	rm -f */*.pyc