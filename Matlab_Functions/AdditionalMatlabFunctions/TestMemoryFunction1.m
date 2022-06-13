function y1 = TestMemoryFunction1(x1)

    x1(1) = 5;
    y1 = TestMemoryFunction2(x1);
    y1(1) = 1;
    y1(2) = 2;

end