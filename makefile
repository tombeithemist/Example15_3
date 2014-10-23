# NB!!! Must be tab characters indeneting the commands, not spaces.

# Under Windows rm does not work when invoked by make, but does work if invoked by nmake.
DEL = rm
SRC = src
INCLUDE = src
OUT = Debug
TARGET = Debug

LINK_TARGET = $(TARGET)/Example15_3.exe

# Here is a Make Macro that uses the backslash to extend to multiple lines.
OBJS =  \
 $(OUT)/Example15_3.o \
 $(OUT)/SparseMatrix.o \
 $(OUT)/Boundary.o 
  
# The rule for "all" is used to incrementally build your system.
# It does this by expressing a dependency on the results of that system,
# which in turn have their own rules and dependencies.
all : $(LINK_TARGET)
	echo All done

# See note above as this only runs under nmake	
clean:
	$(DEL) $(OUT)/Example15_3.o $(TARGET)/Example15_3.exe

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
$(OUT)/%.o : $(SRC)/%.cpp
	g++ -g -o $@ -c $<

# These are Dependency Rules, which are rules without any command.
# Dependency Rules indicate that if any file to the right of the colon changes,
# the target to the left of the colon should be considered out-of-date.
# The commands for making an out-of-date target up-to-date may be found elsewhere
# (in this case, by the Pattern Rule above).
# Dependency Rules are often used to capture header file dependencies.
$(OUT)/Example15_3.o : $(SRC)/Example15_3.cpp $(INCLUDE)/SparseMatrix.h
			