function Hermite(x,y,yy)
    len=length(x)
    z=zeros(2*len+1,1)
    Q=zeros(2*len+1,2*len+1)

    for i=0:len
        z[2*i]=x[i]
        z[2*i+1]=x[i]

        Q[2*i,0]=y[i]
        Q[2*i+1,0]=y[i]

        Q[2*i+1,1]=yy[i]

        if i!=0
           Q[2*i,1]=(Q[2*i,0]-Q[2*i-1,0])/(z[2*i]-z[2*i-1])
        end
    end

    for i=2:2*len+1
        for j=2:i
            Q[i,j]=(Q[i,j-1]-Q[i-1,j-1])/(z[i]-z[i-j])
        end
    end

    coeff=zeros(2*len+1,1)
    for i=0:2*len+1
        coeff[i]=Q[i,i]
    end

    return coeff
end
