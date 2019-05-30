import pytest

def function_to_test(a,b):
    return a+b

def test__function_to_test():
    a = 1
    b = 2

    expected_result = a + b
    assert function_to_test(a,b) == expected_result

def test__function_to_test2():
    a = 1
    b = 2

    expected_result = a + b
    assert function_to_test(a+1,b) == expected_result
