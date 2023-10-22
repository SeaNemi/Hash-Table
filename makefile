CXX = g++
CXXFLAGS = -Wall -g
IODIR =../../proj0_IO/

proj0: vdetect.o mytest.cpp
	$(CXX) $(CXXFLAGS) vdetect.o mytest.cpp -o proj4

vdetect.o: vdetect.h vdetect.cpp
	$(CXX) $(CXXFLAGS) -c vdetect.cpp

clean:
	rm *.o*
	rm *~

run:
	./proj4

val:
	valgrind ./proj4