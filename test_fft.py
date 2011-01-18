a = arange(18).reshape((2,3,3))
aI = rfftn(a)

print rfft(a,axis=1)

b = zeros((2,2,3),'complex')
for i in xrange(2):
    for j in xrange(3):
        b[i,:,j] =  rfft(a[i,:,j])

print b
