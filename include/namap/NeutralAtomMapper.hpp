//
// Created by Ludwig Schmid on 05.10.23.
//

#pragma once

#include "Architecture.hpp"

constexpr std::int16_t DEFAULT_POSITION = -1;

class test {
  int a;

  test(int a) : a(a){};

  int getA() { return a; };

  void setA(int _a) { this->a = _a; };

  void printA() { std::cout << a << std::endl; };
};

int return_1();
