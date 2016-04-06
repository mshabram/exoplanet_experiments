# -*- coding: utf-8 -*-
"""
Created on Fri Apr 01 16:07:42 2016

@author: jcat
"""

def memoize(function):
    cache = {}
    def decorated_function(*args):
        if args in cache:
            return cache[args]
        else:
            val = function(*args)
            cache[args] = val
            return val
    return decorated_function

# show example usage
if __name__ == "__main__":
  import cProfile

  def fib_slow(n):
    "Returns the nth Fibonacci number (slow, non-memoized version)"
    if n == 0: return 0
    if n == 1: return 1
    return fib_slow(n-2) + fib_slow(n-1)

  @memoize
  def fib_fast(n):
    "Returns the nth Fibonacci number (fast, memoized version)"
    if n == 0: return 0
    if n == 1: return 1
    return fib_fast(n-2) + fib_fast(n-1)

  print "We're now going to run two versions of the same function."
  print "Both calculate the 35th Fibonacci number."
  print ""
  print "The first function is not memoized, and thus very slow."
  print "The second is memoized, using our decorator, and thus very fast."
  print ""
  print "Check out the speed differences between the two."
  print ""

  print "Calling fib_slow(35) ..."
  cProfile.run("fib_slow(35)")

  print "Calling fib_fast(35) ..."
  cProfile.run("fib_fast(35)")