class foo(object):
    def __init__(self, x = []):
        if not isinstance(x, list):
            raise TypeError('x must be a list')
        else:
            self.x = list(x)

f1 = foo()
f2 = foo()
f1.x = [1,2]
print 'f1.x', f1.x,' f2.x = ',f2.x

f3 = foo()
f4 = foo()
f3.x.append(2)
print 'f3.x = ',f3.x,' f4.x = ',f4.x
