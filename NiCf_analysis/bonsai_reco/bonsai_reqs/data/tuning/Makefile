ROOTCFLAGS = $(shell root-config --cflags)
ROOTLDFLAGS = $(shell root-config --ldflags)
ROOTLIBS = $(shell root-config --libs)
ROOTINCS = $(shell root-config --incdir)

CXXFLAGS = $(ROOTCFLAGS) -Wall
CXXLIBS = $(ROOTLIBS)
CXXINCS = -I$(ROOTINCS) -I.
LDFLAGS = $(ROOTLDFLAGS) 

all: show_like create_like

binfile.o: binfile.cc binfile.h
	g++ $(CXXFLAGS) -c binfile.cc

create_like: create_like.cc binfile.o
	g++ $(CXXFLAGS) -o create_like -I. create_like.cc binfile.o

show_like: show_like.cc binfile.o
	g++ $(CXXFLAGS) -o show_like $(CXXINCS) show_like.cc binfile.o $(CXXLIBS)
clean:
	rm binfile.o show_like create_like *.bin
