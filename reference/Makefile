# Builds the Kikuchi-Kanamori reference programs with gfortran.
# -std=legacy accepts the Fortran 77 extensions the ETAL code uses.
# The code reads locals before setting them and expects them to keep values
# between calls, as static allocation did on the compilers it was written
# for. -fno-automatic -finit-local-zero restore that; without them plotcm,
# mom3 and gr3 change output from run to run.

FC      ?= gfortran
FFLAGS  ?= -std=legacy -O -fno-automatic -finit-local-zero
BUILD   := build

PROGRAMS := $(basename $(notdir $(wildcard main/*.f)))
SUB_OBJS := $(patsubst subroutine/%.f,$(BUILD)/obj/%.o,$(wildcard subroutine/*.f))
REC_OBJS := $(patsubst recipes/%.f,$(BUILD)/obj/%.o,$(wildcard recipes/*.f))

all: $(addprefix $(BUILD)/bin/,$(PROGRAMS))

$(BUILD)/obj $(BUILD)/bin:
	mkdir -p $@

$(BUILD)/obj/%.o: subroutine/%.f | $(BUILD)/obj
	$(FC) $(FFLAGS) -c $< -o $@

$(BUILD)/obj/%.o: recipes/%.f | $(BUILD)/obj
	$(FC) $(FFLAGS) -c $< -o $@

$(BUILD)/mylib.a: $(SUB_OBJS)
	ar rcs $@ $^

$(BUILD)/recipes.a: $(REC_OBJS)
	ar rcs $@ $^

$(BUILD)/bin/%: main/%.f $(BUILD)/mylib.a $(BUILD)/recipes.a | $(BUILD)/bin
	$(FC) $(FFLAGS) $< $(BUILD)/mylib.a $(BUILD)/recipes.a -o $@

golden: all
	./regenerate-golden

clean:
	rm -rf $(BUILD)

.PHONY: all golden clean
