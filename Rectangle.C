// example: one class, two objects
#include <iostream>
//#include <stdio.h>
using namespace std;

class Rectangle {
    int width, height;
  public:
    void set_values (int,int);
    int area () {return width*height;};
};

void Rectangle::set_values (int x, int y) {
  width = x;
  height = y;
}

int main () {
  Rectangle rect, rectb;
  rect.set_values(3,4);
  rectb.set_values(5,6);
  std::cout << "rect area: " << rect.area() << std::endl;
  std::cout << "rectb area: " << rectb.area() << std::endl;
  return 0;
}
