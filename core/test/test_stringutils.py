from tamasis.stringutils import *

class TestFailure(): pass

if strenum(['blue', 'red', 'yellow'], 'or') != "'blue', 'red' or 'yellow'": raise TestFailure()

if strplural('cat', 0, prepend=False, s='' ) != 'cat' or \
   strplural('cat', 1, prepend=False, s='' ) != 'cat' or \
   strplural('cat', 2, prepend=False, s='' ) != 'cats' or \
   strplural('cat', 0, prepend=True,  s='' ) != 'no cat' or \
   strplural('cat', 1, prepend=True,  s='' ) != '1 cat' or \
   strplural('cat', 2, prepend=True,  s='' ) != '2 cats' or \
   strplural('cat', 0, prepend=False, s=':') != 'cat' or \
   strplural('cat', 1, prepend=False, s=':') != 'cat:' or \
   strplural('cat', 2, prepend=False, s=':') != 'cats:' or \
   strplural('cat', 0, prepend=True,  s=':') != 'no cat' or \
   strplural('cat', 1, prepend=True,  s=':') != '1 cat:' or \
   strplural('cat', 2, prepend=True,  s=':') != '2 cats:': raise TestFailure()
