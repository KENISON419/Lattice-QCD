##################################
TARGET = ising.x

##########

CC  = icc
CXX = icpc

CFLAGS   = -O3
CXXFLAGS = -O3

##########

OBJS = \
	main.o \
	ising.o \
	mc.o \
	measurement.o \
	mt19937ar_lib.o \
	stat_utils.o \

##########

all: $(TARGET)

$(TARGET) : $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $(OBJS)

clean:
	rm -f $(TARGET) ; \
	rm -f *.o ; \

##########

.c.o:
	$(CC) $(CFLAGS) -c -o $@ $<

.C.o:
	$(CXX) $(CXXFLAGS) -c -o $@ $<

#################   EOF   #################
