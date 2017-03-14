CXX = gcc
CXXFLAGS = -fopenmp  -w -O3 
LDFLAGS = -fopenmp -O3 -w
OFLAGS= -O3
WFLAGS= -w
SOURCE = $(guass_p.c)
OBJECTS = $(SOURCE:.c=.o)
TARGET = guass_p


default: $(TARGET)
        

%.o: %.c
	$(CXX)  $(CXXFLAGS) $(OFLAGS) $(WFLAGS) -c -o $@ $< $(LDFLAGS)
        
       

parcount: $(OBJECTS)
	$(CXX) $(CXXFLAGS) $(OFLAGS) $(WFLAGS) $< -o $@ $(LDFLAGS)
       
 
clean:
	rm -f $(OBJECTS) $(TARGET)
        


