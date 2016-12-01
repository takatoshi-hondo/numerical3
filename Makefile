CXX=g++
AR=ar
CXXFLAGS+=-L/usr/local/lib -I/usr/local/include -g
LDLIBS=-lm
OBJS=n3matrix.o n3complexnumber.o
TARGET=libn3.a
EXE=demo

.cpp.o:
	$(CXX) $(CXXFLAGS) -c $^ -o $@ $(LDLIBS)

all: $(OBJS) $(TARGET) 

$(TARGET): $(OBJS)
	$(AR) rv $(TARGET) $(OBJS)

$(EXE): $(EXE).cpp $(OBJS)
	$(CXX) $(CXXFLAGS) $(EXE).cpp -o $@ -ln3

install:
	@install -C $(TARGET) /usr/local/lib
	mkdir -p /usr/local/include/numerical3
	@install -C *.hpp /usr/local/include/numerical3

clean:
	rm -f $(TARGET) $(OBJS) $(EXE) *~

