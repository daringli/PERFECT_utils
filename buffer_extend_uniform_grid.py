import numpy

def buffer_extend_uniform_grid(xmin,xmax,dx,first=True):
    Nx = int(numpy.floor((xmax - xmin)/dx))
    x = numpy.array([xmin + i*dx for i in range(Nx+1)])
    if not first:
        x = x[1:]
    return x

if __name__=="__main__":
    #check that it works as expected
    
    xmin = 0
    xmax = 3.5764
    dx = 1.0
    ans = [0,1,2,3]
    assert((buffer_extend_uniform_grid(xmin,xmax,dx) == numpy.array(ans)).all())

    xmin = 0
    xmax = 4.0
    dx = 1.0
    ans = [0,1,2,3,4]
    assert((buffer_extend_uniform_grid(xmin,xmax,dx) == numpy.array(ans)).all())

    xmin = 1.0
    xmax = 4.0
    dx = 1.0
    ans = [1,2,3,4]
    assert((buffer_extend_uniform_grid(xmin,xmax,dx) == numpy.array(ans)).all())

    xmin = 0.645
    xmax = 4.0
    dx = 1.0
    ans = [0.645,1+0.645,2+0.645,3+0.645]
    assert((buffer_extend_uniform_grid(xmin,xmax,dx) == numpy.array(ans)).all())

    xmin = 0.645
    xmax = 3.25
    dx = 1.0
    ans = [0.645,1+0.645,2+0.645]
    assert((buffer_extend_uniform_grid(xmin,xmax,dx) == numpy.array(ans)).all())
    
    xmin = 0.645
    xmax = 3.25
    dx = 1.25
    ans = [0.645,1.25+0.645,1.25+1.25+0.645]
    assert((buffer_extend_uniform_grid(xmin,xmax,dx) == numpy.array(ans)).all())
    
