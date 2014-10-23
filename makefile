# NB!!! Must be tab characters indeneting the commands, not spaces.

# Under Windows rm does not work when invoked by make, but does work if invoked by nmake.
DEL_CMD = rm
SRC_DIR = src
INCLUDE_DIR = src
OUT_DIR = Debug
TARGET_DIR = Debug

LINK_TARGET = $(TARGET_DIR)/Example15_3.exe

# Here is a Make Macro that uses the backslash to extend to multiple lines.
OBJS =  $(OUT_DIR)/Boundary.o \
        $(OUT_DIR)/Example15_3.o \
        $(OUT_DIR)/SparseMatrix.o
  
# The rule for "all" is used to incrementally build your system.
# It does this by expressing a dependency on the results of that system,
# which in turn have their own rules and dependencies.
all : $(LINK_TARGET)
	echo All done

# See note above as this only runs under nmake	
clean:
	$(DEL_CMD) $(OBJS) $(LINK_TARGET)

# Here is a Rule that uses some built-in Make Macros in its command:
# $@ expands to the rule's target.
# $^ expands to the rule's dependencies
$(LINK_TARGET) : $(OBJS)
	g++ -g -o $@ $^

# Here is a Pattern Rule, often used for compile-line.
# It says how to create a file with a .o suffix, given a file with a .cpp suffix.
# The rule's command uses some built-in Make Macros:
# $@ for the pattern-matched target
# $< for the pattern-matched dependency
$(OUT_DIR)/%.o : $(SRC_DIR)/%.cpp
	g++ -g -o $@ -c $<

# These are Dependency Rules, which are rules without any command.
# Dependency Rules indicate that if any file to the right of the colon changes,
# the target to the left of the colon should be considered out-of-date.
# The commands for making an out-of-date target up-to-date may be found elsewhere
# (in this case, by the Pattern Rule above).
# Dependency Rules are often used to capture header file dependencies.
# This does not seem to work as expected??
# $(OUT_DIR)/Example15_3.o : $(INCLUDE_DIR)/SparseMatrix.h $(INCLUDE_DIR)/Boundary.h
			