// log_heap.cpp
#include <malloc.h>
#include <iostream>

void log_heap() {
    auto info = mallinfo2();
    std::cout << "Heap in use: " << info.uordblks / 1024 << " KB\n";
}
