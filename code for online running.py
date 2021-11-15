from matplotlib import pyplot as plt
import math as m

#________________________GIVEN DATA________________________________________

print("-------------------INPUT DATA-------------------- \n")

print("-----INPUT DATA FOR CHANNEL GEOMETRY-----\n")
try:
    Q=float(input("Enter the value for discharge,in in m^3/s:")) #discharge      
    B=float(input("Enter the value for channel width,in in m:")) #bottom width       
    Z1=float(input("Enter the value for side slope-1(m),(m:1):")) #side slope 1
    Z2=float(input("Enter the value for side slope-2(n),(n:1):")) #side slope 2
    n=float(input("Enter the value for Manning's coefficient:")) #Mannings coefficient
    So=float(input("Enter the value for bed slope: ")) #bed slope

    print("\n-----INPUT DATA FOR INITIAL CONDITIONS OF THE FLOW-----\n")

    x0=float(input("Enter the initial value for x,in metres(for GVF calculation):"))    #initial condition
    y0=float(input("Enter the initial value for y,in metres(for GVF calculation):"))    #initial condition    
    Y1=y0 #known depth    in m
    Xm=float(input("Enter the maximum value for x over which the flow is to be observed,in metres(for GVF calculation):")) #maximum distance   in m
    distanceX=float(input("Enter the distance x where the flow depth is to be calculated,in metres(for GVF calculation):")) #desired distance distance   in m
    distanceX=int(distanceX)

    print("\n--------------------OUTPUT DATA-------------------- \n")

#_____________________CALCULATION OF NORMAL DEPTH__________________________

# Newton Raphson method 
    def func(Y): 
        return ((2*B+(Z1+Z2)*Y)*Y*0.5)**5-(Q*n/m.sqrt(So))**3*((B+(m.sqrt(Z1**2)+m.sqrt(Z2**2))*Y))**2
    
    def derivFunc( Y ): 
        return 5*(( 2*B+(Z1+Z2)*Y)*Y*0.5)**4*(B+Y*(Z1+Z2))-(Q*n/m.sqrt(So))**3*2*(B+(m.sqrt(Z1**2)+m.sqrt(Z2**2))*Y)*(m.sqrt(Z1**2)+m.sqrt(Z2**2))

    def newtonRaphson( Y ): 
        h = func(Y) / derivFunc(Y) 
        while abs(h) >= 0.0001: 
            h = func(Y)/derivFunc(Y)        
            Y= Y - h     
        print("The NORMAL depth is(in m) : ", "%.4f"% Y)
        return Y
    Yn=newtonRaphson(Y1)

#________________________CALCULATION OF CRITICAL DEPTH___________________________

# Newton Raphson method 
    def func2( Y): 
        return ((2*B+(Z1+Z2)*Y)*Y*0.5)**3-(Q**2/9.8)*((B+(Z1+Z2)*Y))

    def derivFunc2( Y ): 
        return 3*((2*B+(Z1+Z2)*Y)*Y*0.5)**2*(2*B+2*(Z1+Z2)*Y)-(Q**2/9.8)*(Z1+Z2) 
 
    def newtonRaphson2( Y ): 
        h = func2(Y) / derivFunc2(Y) 
        while abs(h) >= 0.0001: 
            h = func2(Y)/derivFunc2(Y)          
            Y= Y-h 
        print("The CRITICAL depth is (in m) : ",  "%.4f"% Y)
        return Y
    Yc=newtonRaphson2(Y1)


#________________________PLOTTING OF THE BED SLOPE______________________________

#uu and vv are empty lists that stores the co ordinates  of the bed slope
    uu=[]
    vv=[]
    h=1   #step size
    z=int(Xm/h+1)    #no of steps
    u=0
    for i in range(z):
        uu.append(u)
        v=So*(Xm-u)
        vv.append(v)
        u=u+h

    plt.plot(uu,vv,label='CHANNEL BED')

#_____________________________PLOTTING OF NDL AND CDL__________________________
#r stores the co ordinates of ndl....R is an empty list
#s stores the co ordinates of cdl.....S is an  empty list
    R=[ ]
    S=[ ]
    for i in range(z):
        r=vv[i]+Yn
        R.append(r)
        s=vv[i]+Yc
        S.append(s)
    plt.plot(uu,R,label='NDL')
    plt.plot(uu,S,label='CDL')

#__________________GRADUALLY VARIED FLOW PARAMETERS__________________________
#GIVEN

    h=1     #step size
    a=[]    #To store x-coordinate of GVF flow
    b=[]    #To store y-coordinate of GVF flow+y-coordinate of bed slope
    c=[]    #To store x-coordinate of GVF flow
    a.append(x0)
    b.append(y0)

#________________________________EULER'S METHOD______________________


    def func( x, y ):
        t=B+(Z1+Z2)*y
        ar=(B+t)*y*0.5
        p=B+(m.sqrt(Z1**2+1)+m.sqrt(Z2**2+1))*y
        f1=4/3
        f2=10/3
        Sf=((n**2)*(Q**2)*(p**f1)/(ar**f2))
        return((So-Sf)/(1-(t*Q**2)/(9.8*ar**3)))

# Function for euler formula

    def euler( x0, y, h, x1):
        j=0
        while x0 < x1:       
            y = y + h * func(x0, y)
            x0 = x0 + h
            a.append(x0)
            c.append(y)          # y co-ordinate excluding bed slope
            b.append(vv[j]+y)    # y co-ordinate including bed slope
            j=j+1
            if distanceX==x0:
                print("Approximate solution at x(in m) = ", x0, " is ", "%.6f"% y)    # Printing approximate value at desired distance
     
    euler(x0, y0, h, Xm)

#_______________________PLOTTING OF GVF PROFILE_________________

#_____DETERMINING THE FLOW PROFILE____________

    if (Yn>Yc):                                                 # for MILD Slope
        if (c[10]>Yn and c[10]>Yc):
            print('\n Type of GVF flow profile is:  M1 \n')
        elif (c[10]<Yn and c[10]>Yc):
            print('\n Type of GVF flow profile is:  M2 \n')
        else:
            print('\n Type of GVF flow profile is:  M3 \n')
    if (Yn<Yc):                                                 # for STEEP Slope
        if(c[10]>Yn and c[10]>Yc):
            print('\n Type of GVF flow profile is:  S1 \n')
        elif(c[10]>Yn and c[10]<Yc):
            print('\n Type of GVF flow profile is:  S2 \n')
        else:
            print('\n Type of GVF flow profile is:  S3 \n')
    if (Yn==Yc):                                                # for CRITICAL Slope
        if(c[10]>Yn and c[10]>Yc):
            print('\n Type of GVF flow profile is:  C1 \n')
        else:
            print('\n Type of GVF flow profile is:  C3 \n')
    uu.pop(0)
    b.pop(0)
    plt.plot(uu,b,label='WATER SURFACE')
    plt.xlabel('x in metres')
    plt.ylabel('y in meres')
    plt.title('Profile of gradually varied flow')
    plt.legend()
    plt.show()
except:
    pass


    

