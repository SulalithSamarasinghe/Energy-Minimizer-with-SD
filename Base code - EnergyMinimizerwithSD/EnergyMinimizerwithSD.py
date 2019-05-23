#
#Energy Minimizer with SD is a programme that use Steepest Desents method to
#find the minimum value of a given function.
#
import sympy as sym
import csv
#
x,y,z,c,dx,dy = sym.symbols('x y z c dx dy')
#
print('+++++++++++++++++++++++++')
print('     Energy Minimizer    ')
print('           with          ')
print(' Steepest Desents Method ')
print('+++++++++++++++++++++++++')
print('')
print('Find the minimum value of a function with Energy Minimizer!')
print('')
#
sf = input('Enter the function to be minimized : ')
sf = eval(sf)
gf = (sf.diff(x))+(sf.diff(y))
gfX = sf.diff(x)
gfY = sf.diff(y)
#
print('')
spX = float(input('Enter the starting point x coordinate : '))
spXp = spX
print('')
spY = float(input('Enter the starting point y coordinate : '))
spYp = spY
print('')
fileName = input('Enter the output file name : ')
#
with open(fileName,'w',newline='') as dFile:
    thewriter = csv.writer(dFile)
    thewriter.writerow([])
    thewriter.writerow(['Function used to minimize = ',sf])
    thewriter.writerow(['Gradient of the Function  = ',gf])
    thewriter.writerow(['First derivative of the function (df1)  = ',(sf.diff(x)*dx)+(sf.diff(y)*dy)])
    thewriter.writerow([])
    thewriter.writerow(['Iteration','Starting point x coordinate','Starting point y coordinate','Direction of the force x coordinarte','Direction of the force y coordinarte','Gradient of the line','Equation of the line (f2)','df2','df1 + z*df2','Value of z','New x coordinate','New y coordinate','Minimum point value'])
#
iter = 1
#
def SteepestDesentsCalculation(sf,gfX,gfY,spX,spY,spXp,spYp,iter):
    #
    pmV = sf.subs(x,spX)
    pmV = pmV.subs(y,spY)
    pmV = float(pmV)
    #
    c = sym.symbols('c')
    #
    gX = gfX.subs(x,spX)
    gY = gfY.subs(y,spY)
    #
    dfX = (-1)*gX
    dfY = (-1)*gY
    #
    gl = dfY/dfX
    #
    el = x + c - (y*(1/gl))
    #
    elS = el.subs(x,spX)
    elS = elS.subs(y,spY)
    #
    c = sym.solve((elS),(c))
    c = float(c[0])
    el = x + c - (y*(1/gl))
    #
    df1 = (sf.diff(x)*dx)+(sf.diff(y)*dy)
    df2 = (el.diff(x)*dx)+(el.diff(y)*dy)
    #
    lum = df1 + z*df2
    lump = df1 + z*df2
    lum = sym.expand(lum)
    #
    dfx = lum.subs(dy,0)
    dfx = dfx.subs(dx,1)
    #
    spX = sym.solve((dfx),(x))
    spX = spX[0]
    spXf = x + spX
    #
    dfy = lum.subs(dx,0)
    dfy = dfy.subs(dy,1)
    #
    spY = sym.solve((dfy),(y))
    spY = spY[0]
    spYf = y + spY
    #
    elS = el.subs(x,spX)
    elS = elS.subs(y,spY)
    #
    zv = sym.solve((elS),(z))
    zv = zv[0]
    #
    spX = float(spX.subs(z,zv))
    spY = float(spY.subs(z,zv))
    #
    fmV = sf.subs(x,spX)
    fmV = fmV.subs(y,spY)
    fmV = float(fmV)
    with open(fileName,'a',newline='') as dFile:
        thewriter = csv.writer(dFile)
        thewriter.writerow([iter,spXp,spYp,dfX,dfY,gl,el,df2,lump,zv,spX,spY,fmV])
    difference = pmV - fmV
    if(difference == 0):
        print('')
        print('Your calculation is over!!!')
        print('')
    else:
        iter += 1
        SteepestDesentsCalculation(sf,gfX,gfY,spX,spY,spX,spY,iter)
    #
SteepestDesentsCalculation(sf,gfX,gfY,spX,spY,spXp,spYp,iter)
