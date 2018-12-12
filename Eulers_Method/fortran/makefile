#Nina Haas
#Makefile for main.F90 and Euler.F90

OBJECTS = main.o euler.o
.PHONY: clean help

main.exe: $(OBJECTS)
	gfortran $(OBJECTS) -o  main.exe

%.o : %.F90
	gfortran -c -g -fdefault-real-8 $<

clean:
	rm -f $(OBJECTS) main.exe

help:
	@echo "Valid targets:"
	@echo "	main.exe"
	@echo "	main.o"
	@echo "	euler.o"
	@echo "	clean: removes .o and .exe files"
