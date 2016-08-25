PROG := chev
SRCS :=	chev.cpp main.cpp

OBJS :=	$(SRCS:%.cpp=%.o)
DEPS :=	$(SRCS:%.cpp=%.d)

VPATH  = src src/chebyshevpolynomialbdg
CXX = g++
CXXFLAGS = -Wextra -O3 -pipe -std=c++14 -Isrc/eigen-eigen-dc6cfdf9bcec

all: $(PROG) ;
#rm -f $(OBJS) $(DEPS)

-include $(DEPS)

$(PROG): $(OBJS)
		$(CXX) $(LDFLAGS) $(CXXFLAGS) -o $@ $^

%.o: %.cpp
		$(CXX) $(CXXFLAGS) -c -MMD -MP $<

clean:
		rm -f $(PROG) $(OBJS) $(DEPS)
