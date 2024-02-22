#include "Point.h"
#include <iostream>
#include <vector>

int main(int argc, char *argv[])
{
    Point pt(5);
    std::cout << pt << std::endl;

    Point pt2({1,2,3});

    std::cout << pt2 << std::endl;

    std::vector<Point> points = Point::random_points(10, 3, 134);

    for (auto itr = points.begin(); itr != points.end(); itr++)
        std::cout << *itr << std::endl;

    std::vector<Point> points2 = Point::random_points(20, 5);

    for (auto itr = points2.begin(); itr != points2.end(); itr++)
        std::cout << *itr << std::endl;

    return 0;
}
