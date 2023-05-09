#---------------- Lagrange's Polynomial------------------------------
#--------------------------------------------------------
x = [1 2 4 5 7] # iteration through this x value
y1 = [52 5 -5 -40 10]
x_v = 3 # f(xv) prediction
y = 0;
l = 1;
function langrage(x,y1,x_v,y)
for i in 1:length(x)

for j in 1:length(x)
    if(j!=i)
    l = l*((x_v-x[j])*(x[i]-x[j]).^-1)
    end
end
y = y+(y1[i].*l)
display(y)
l = 1;
end
end


#------------------------------------Langrage poly--------------------------------------
x = [1 2 4 5 7] # iteration through this x value
y1 = [52 5 -5 -40 10]
x_v = 3 # f(xv) prediction
y = 0;
l = 1;
using SymPy
b = symbols("b")
a = symbols("a")
function langrage1(x,y1,y,l)
for i in 1:length(x)

for j in 1:length(x)
    if(j!=i)
    l = l*((a-x[j])*(x[i]-x[j]).^-1)
    end
end
display(l)

y = y+(y1[i].*l)
l = 1;

end
y = simplify(y)
display(y)
display(" ")
end

#----------------------------Divided Differences--------------------------
x = [0.0 20.0 40.0 60.0 80.0 100.0]
y1 = [26.0 48.6 61.6 71.2 74.8 75.2]
#----------------
#example ->
x = [0 2 3]
y = [1 2 4]
#--------------
N1 = 0;
N2 = 0;
N3 = 0
N4 = 0;
a = 0;
a2= 0;
a3 = 0;
xv = 55;


for i in 1:(length(x)-1)
        a = (y1[i+1]-y1[i])*(x[i+1]-x[i]).^-1
        display(a)
        N1 = [N1 a]
end
for j in 1:(length(N1)-2)
a2 = (N1[j+2]-N1[j+1])*(x[j+2]-x[j]).^-1
display(a2)
N2 = [N2 a2]
end

for m in 1:(length(N2)-2)
a3 = (N2[m+2]-N2[m+1])*(x[m+3]-x[m]).^-1
display(a3)
N3 = [N3 a3]
end

for g in 1:(length(N3)-2)
a4 = (N3[g+2]-N3[g+1])*(x[g+4]-x[g]).^-1
display(a4)
N4 = [N4 a4]
end
N5 = 0
a5 = 0
for g in 1:(length(N4)-2)
a5 = (N4[g+2]-N4[g+1])*(x[g+5]-x[g]).^-1
display(a5)
N5 = [N5 a5]
end
Nf = y1[1]
Nf = [Nf N1[2] N2[2] N3[2] N4[2] N5[2]]
l = 1;
s = 1
for i in 1: 5
 l = l*(xv-x[i])
 s = [s l]
end
P = 0;
for i in 1 : 6
 P = P+Nf[i]*s[i]
end






#----------------------------Newton Method 2nd ----------------------------------------------
x = [0 1 2 3]
y = [2 1  0 -1 ]

#---------------------------------------------

using SymPy

x = [-0.1 0.0 0.2 0.3 0.35 0.55]
y = [5.30000 2.00000 3.19000 1.00000 0.97260 0.23322]

                    #
#------------------------------Newton Divided Diff.-----------------#
function NewtonDividedDiff(x,y)

y1 = y
F = y1[1];
xv = 0;
N = 0;
l = length(x)
m = length(x)-1;
h = 1
keep_h = 1;

for i in 1:l
    print(y1[i])
    print(" || ")
end
display("->")
for i in 1 : m

    for i in 1:(l-1)
        N1 = round((y1[i+1]-y1[i])*(x[h+1]-x[i]).^-1,digits=10)
        N = [N N1]
        print(N1)
        print(" || ")

        h = h+1;
    end
display("->")

    h = keep_h
    keep_h = keep_h+1;
    h = h+1;

    y1 = N[2]
    for j in 3:(length(N))
     y1 = [y1 N[j]]
    end


    F = [F y1[1]]
    N = 0 #resetting N
    l = l-1;

end

a = symbols("a")
val = 1
P1 = 0

for i in 1:(length(F))
    P1 = P1 + F[i]*val
    val = val*(a-x[i])


end
println(P1)
P1 = simplify(P1)
println(P1)
y_output(a) = P1(a);
y_outputPrint = y_output.(x)
return y_outputPrint
end

y_P = NewtonDividedDiff(x,y)
using Plots
plotly()
scatter(x,y)
x1= vec(x)
y_P = vec(y_P)
plot!(x1,y_P, linewidth=3)

#additional Part -> plotting from obtained poly
f1(a) = -7318.2347134348*a^5 + 8218.9194213192*a^4 - 2746.49922817896*a^3 + 194.506847078794*a^2 + 22.8664198823317*a + 2.0
y1 = f1.(x)
x = vec(x)
y1 = vec(y1)
plot!(x,y1,lw=3)

# -------------Calculating or rather predicting the value --------------------

b = 1;
s = 1
for i in 1: length(y)-1
 b= b*(xv-x[i])
 s = [s b]
end
P = 0;
for i in 1 : length(y)
 P = P+F[i]*s[i]
end
#-------------------------------Chebyshev's nodes--------------------------------------------
n = 3; #deg of the polynomial
x = 0;
k = 1;
b = 0;

function ChebyshevsNodes(n,k,a,b)

k = 1;
x = 0;
m = n+1;

for i in 1:m
    l = cos.(((2*k-1)*pi)*(2*n).^-1)
    x = [x l]
    k = k+1;
end
display(x)
x_new=0;
g = 0;
x1 = x[2]
for h in 3:length(x)
    x1 = [x1 x[h]]
end
display(x1)
for j in 1:length(x1)

    x_new = 1/2*((b-a)*x1[j]+a+b)
    g = [g x_new]
end
display(g)
p = g[2]
    for d in 3:length(g)
        p = [p g[d]]
    end
display(p)
return p
end



#-------------------------------------Chebyshevs polynomial re version-----------------#


function ChebyshevsNodes(n,a,b,f)

k = 1;
x = 0;
m = n+1;


for i in 1:m
    l = cos.(((2*k-1)*pi)*(2*n).^-1)
    x = [x l]
    k = k+1;
end
display(x)
x_new=0;
g = 0;
x1 = x[2]
for h in 3:length(x)
    x1 = [x1 x[h]]
end
display(x1)
for j in 1:length(x1)

    x_new = 1/2*((b-a)*x1[j]+a+b)
    g = [g x_new]
end
display(g)
p = g[2]
    for d in 3:length(g)-1
        p = [p g[d]]
    end
display(p)


y1 = f.(p)
y111 = NewtonDividedDiff(p,y1)

return p,y111
end


#------------------Example -----
#f(x) = x*e^(x) [-5,5]

x = -5:0.1:5
y_fun_x = x.*exp.(x)


#----------------Using Cheby approx-----------#
f11(x1) = x1.*exp.(x1)
x1 = ChebyshevsNodes(15,-5,5,f11)[1]
y1 = ChebyshevsNodes(15,-5,5,f11)[2]


using Plots
plotly()

plot(x,y_fun_x,linewidth=3)
scatter!(x1,y1)






#----------------------------------------------

p1 = ChebyshevsNodes(2,-1,1)

x = [ 0.933013  0.5  0.0669873 ]
y= 1/2*cos.(x)+1/3*sin.(2*x)

langrage1(x,y,0,1)

NewtonDividedDiff(x,y)

#answer
# -0.742773223*a^2 + 0.827316500510899*a + 0.491316664708963
#---------------------------
f = x-> exp.((sin.(x)).^3)+x.^6+2*x.^4-x.^3-1

NewtonMethod(f,derivative(f),-2)
BisectionMethod(f,-2,2)
f = x->1/2+(1/4)*x^2-x*sin.(x)-1/2*cos.(2*x)
#--------------------------------


#----------------------------------------

x = [ 0.92388  0.382683  -0.382683  -0.92388 ]
x = [-5 -4 -3 -2 -1 0 1 2 3 4 5]
y = [5 5 5 5 5 5 5 5 5 5 45]

using PolyFit
fit(x,y)
x = vec(x)
y =sin.(x)
plot(x,y,linewidth= 3)
g=a->-0.1585048895*a^3 + 1.12600206936264e-11*a^2 + 0.998982837200365*a + 1.31725880203604e-12
#error = (f(x)-p(x)) = 0.0002527688...
#----------------

HermiteInterpolationWithDD(x,y,der)









#--------------------------------------
x = 0:0.5:1.5
y = x.*exp.(x)
y = log.(x.+2)
langrage1(x,y,0,1)
NewtonDissection(z,y)
p22 = ChebyshevsNodes(11,1,4,0)
y = z.^-1
z = [ 2.86603  2.0  1.13397 ]
y = (z.+2).^-1
x_new = [1.44291  1.03701  0.462987  0.0570904 ]
x = [0 ,(3)^(1/2)/2 ,-(3)^(1/2)/2]
HermiteInterpolationWithDD(x,y,der)

#-------------------------------------
x = [5.0 ,4.0 ,3.0, 2.0, 1.0, 0 ,-1.0, -2.0, -3.0, -4.0, -5.0]
y  =[45.0 ,5.0, 5.0,5.0, 5.0, 5.0 ,5.0 ,5.0 ,5.0 ,5.0 ,5.0]

using Polynomials
fit(x,y)
NewtonDissection(x,y)
#------------------------------------
g = a->1.10229e-5*a^10 + 5.51148e-5*a^9 - 0.0003306872*a^8 - 0.0016534433*a^7 + 0.00300925549999999*a^6 + 0.0150463255*a^5 - 0.00903878059999932*a^4 - 0.045194108300007*a^3 +4 - 0.045194108300007*a^3 + 0.00634921540000346*a^2 + 0.0317460793000434*a + 5.00000000599998;


#----------------------------

#----------------------Hermite Interpolation requires revision--------------------------
using SymPy
x =  [-1 0 1 ];
y = [2 1 2]
der = [-8 0 8];

function HermiteInterpolationWithDD(x,y,der)
    l = 0;
    l1 = 0;

    for i in 1:length(x)

        for j in 1:2
          l = [l x[i]]
          l1 = [l1 y[i]]
        end

    end
    z = l[2]
    fz = l1[2]

    for i in 3:length(l)
    z = [z l[i]]
    fz  = [fz l1[i]]

    end

    g = 1;
    h = 1;
    N = 0;
    N1 = 0;
    keep_h = 2;
    number = 3;
    F = fz[1];

    for m in 1:length(z)-1
        for i in 1:length(fz)-1

            if(fz[i+1]!=fz[i])
                N1 = round((fz[i+1]-fz[i])*(z[h+1]-z[i]).^-1,digits=10)
                N = [N N1]
                print(N1)
                print(" || ")
            end
            if(fz[i+1]==fz[i])
                N1 = der[g]
                N = [N N1]
                print(N1)
                print(" || ")

                g = g+1;
            end
            h = h+1;
          end
        h = keep_h
        keep_h = keep_h+1;
        fz = N[2]

        for i in 3:length(N)
            fz = [fz N[i]]
        end
        F = [F fz[1]]
        N = 0;
    end
    a = symbols("a")
    val = 1
    P1 = 0

    for i in 1:(length(F))
        P1 = P1 + F[i]*val
        val = val*(a-z[i])
    end
    println(P1)
    P1 = simplify(P1)
    println(P1)
end

HermiteInterpolationWithDD(x,y,der)
#------------------------------Hermite Interpolation with Div Diff-----------------------------------#
using SymPy


function Hermit(x,y,y_p)

    z = x[1]
    j = 1;
    for i = 2:(2*length(x))

        z = [z x[j]]
        if(i%2==0)
            j = j+1
        end

    end

    fz = y[1]
    j = 1;
    for i = 2:(2*length(y))

        fz = [fz y[j]]
        if(i%2==0)
            j = j+1
        end

    end


    F = y_p[1];
    F1 = F;
    Fnew = 0;
    x_t = 2;
    y_t = 2;
    y_pt = 2;

    kx_t = 2;
    ky_t = 2;
    ky_pt = 2;
    n = 1
    i = 2;  #increment for x_t
    p_t = 2
    print(F[1], " || ")

    for m = 1:(length(z)-1)

        for m1 = n:(length(z)-1)

            if(m==1 && m1%2==0)
                F = [F (fz[y_t+1]-fz[y_t])/(z[x_t+1]-z[x_t])]
                print(F[m1], " || ")
                x_t = x_t+2;
                y_t = y_t+2;

            end

            if(m==1 && m1%2!=0 && y_pt<=length(x) && m1!=1)
                F = [F y_p[y_pt]]
                y_pt = y_pt+1
                print(F[m1], " || ")

            end
            if(m!=1)
                Fnew = [Fnew round((F[y_t+1]-F[y_t])/(z[x_t+i]-z[x_t]),digits=10)]
                x_t = x_t+1;
                y_t = y_t+1;


            end

        end


        if(m>=2)
            F1 = [F1 Fnew[2]]
            F = Fnew;
            i = i+1
            display("->")
            F = Fnew[2]
            for j = 2:length(Fnew)
                print(Fnew[j], " || ")
                if(j!=length(Fnew))
                    F = [F Fnew[j+1]]
                end
            end
            Fnew = 0;
        end

        p_t  = 2
        n = n+1;
        x_t = 1
        y_t = 1;

    end


    a = symbols("a")
    val = 1
    P1 = 0

    F1 = [fz[1] F1]

    for i in 1:(length(F1))
        P1 = P1 + F1[i]*val
        val = val*(a-z[i])

    end
    display("----------------------")
    println(P1)
    P1 = simplify(P1)
    println(P1)

    f3 = P1.(z)
    return f3,z

end
#---------------------------Example----------------------#
x =  [1 3 4 10.5];
y = [10 20 3.2 3.552]
y_p = [-3.236256 -3.562 23.1 50];

f111 = Hermit(x,y,y_p)[1]
z1 = Hermit(x,y,y_p)[2]

plot!(vec(z1),vec(f111))


plot(vec(x),vec(y),linewidth=3)
f11(a) =0.1426607989*a^7 - 4.71708327435*a^6 + 59.6046548405*a^5 - 367.9212055517*a^4 + 1199.4328555971*a^3 - 2060.64133849795*a^2 + 1720.7132762627*a - 536.6138201752
scatter!(x,f11.(x))


#-------------------------------


y(x) = (5/(2.75-2.5))*(x-2.5)





















g = a->-0.0027746905*a^5 + 0.02403178315*a^4 - 0.014556057735*a^3 - 0.2352162058455*a^2 - 0.00822919472220004*a + 1.001944055605;
f = x-> 14*x*exp.(x-2)-12*exp.(x-2)-7*x.^3+20x.^2-26*x+12
NewtonMethod(f,derivative(f),0)

f1 = derivative(f)
f2 = derivative(f1)



#-----------------------------------------------------------------------------------------------
#-------------------------Real numbers storing in Julia------------------------------------

el = bitstring(10.375)
el2 = bitstring(-3.465)
el3= bitstring(0.0)
el4 = bitstring(-0.0)

#--------------------------Storing integer values -----------------------------
bitstring(2)
bitstring(-2)

bitstring(0.00) == bitstring(0.0)
bitstring(0)
bitstring(0.0)


#_--------------------------------------------------------------------------------------------

using Interact
using Plots,Colors
timer = Observable(0.0)
@async while true
    sleep[0.1]
    timer[]=timer[]+0.1
end

@manipulate for col1 ="red",col2 = "blue", lw = 1:0.1:10,t = timer
    x = -pi:0.1:pi
    plot(x,[sin.(x.+t).*cos.(x.+2*t)])
end

ui = button()
display(ui)

using Blink
w = Window()
body!(w, "Hello World")
@js w Math.log(10)

@js w console.log("hello, web-scale world")

@js w console.log("hello, web-scale world")

w = Window()
body!(w, ui);


o = observe(ui)

o[]
on(println, o)

filepicker() |> display # observable is the path of selected file
textbox("Write here") |> display # observable is the text typed in by the user
autocomplete(["Mary", "Jane", "Jack"]) |> display # as above, but you can autocomplete words
checkbox(label = "Check me!") |> display # observable is a boolean describing whether it's ticked
toggle(label = "I have read and agreed") |> display # same as a checkbox but styled differently
slider(1:100, label = "To what extent?", value = 33) |> display # Observable is the number selected

dropdown(["a", "b", "c"]) |> display # Observable is option selected
togglebuttons(["a", "b", "c"]) |> display # Observable is option selected
radiobuttons(["a", "b", "c"]) |> display # Observable is option selected
import Colors
@widget wdg function mycolorpicker()
    :r = slider(0:255, label = "red")
    :g = slider(0:255, label = "green")
    :b = slider(0:255, label = "blue")
    @output! wdg Colors.RGB($(:r)/255, $(:g)/255, $(:b)/255)
    @display! wdg plot(sin, color = $(_.output)) ## choose how to display the output, optional
    @layout! wdg hbox(_.display, vbox(:r, :g, :b)) ## custom layout: by default things are stacked vertically
end

mycolorpicker()

using Interact

width, height = 700, 300
colors = ["black", "gray", "silver", "maroon", "red", "olive", "yellow", "green", "lime", "teal", "aqua", "navy", "blue", "purple", "fuchsia"]
color(i) = colors[i%length(colors)+1]
ui = @manipulate for nsamples in 1:200,
        sample_step in slider(0.01:0.01:1.0, value=0.1, label="sample step"),
        phase in slider(0:0.1:2pi, value=0.0, label="phase"),
        radii in 0.1:0.1:60
    cxs_unscaled = [i*sample_step + phase for i in 1:nsamples]
    cys = sin.(cxs_unscaled) .* height/3 .+ height/2
    cxs = cxs_unscaled .* width/4pi
    dom"svg:svg[width=$width, height=$height]"(
        (dom"svg:circle[cx=$(cxs[i]), cy=$(cys[i]), r=$radii, fill=$(color(i))]"()
            for i in 1:nsamples)...
    )
end


x = y = 0:0.1:30

freqs = OrderedDict(zip(["pi/4", "π/2", "3π/4", "π"], [π/4, π/2, 3π/4, π]))

mp = @manipulate for freq1 in freqs, freq2 in slider(0.01:0.1:4π; label="freq2")
    y = @. sin(freq1*x) * sin(freq2*x)
    plot(x, y)
end

using Pluto
Pluto.run()

cpt = OrderedDict(:vertical => Observable(true), :b => slider(1:100), :c => button());
@layout! t Widgets.@map &(:vertical) ? vbox(:b, CSSUtil.vskip(1em), :c) : hbox(:b, CSSUtil.hskip(1em), :c);

a = 30

d = OrderedDict(:label => "My label", :button => button("My button"))
w = Widget{:mywidget}(d)

d = OrderedDict(:label => "My label", :button => button("My button"))
w = Widget{:mywidget}(d)

d = OrderedDict(:label => "My label", :button => button("My button"))
output = map(t -> t > 5 ? "You pressed me many times" : "You didn't press me enough", d[:button])
w = Interact.Widget{:mywidget}(d, output = output)

@layout! w hbox(vbox(:label, :button), observe(_)) # observe(_) refers to the output of the widget

t = Widget{:test}(OrderedDict(:b => slider(1:100), :c => button()));
 @layout! t hbox(:b, CSSUtil.hskip(1em), :c);
@layout! t hbox(:b, CSSUtil.hskip(1em), :c);


cpt = OrderedDict(:vertical => Observable(true), :b => slider(1:100), :c => button());

@js_ win for i in 1:x console.log(i) end


Window()
Window(electron_options::Dict; async=true)

@js win expr

function mycolorpicker()
    r = slider(0:255, label = "red")
    g = slider(0:255, label = "green")
    b = slider(0:255, label = "blue")
    update = button("Update plot")
    output = Interact.@map (&update; Colors.RGB(r[]/255, g[]/255, b[]/255))
    plt = Interact.@map plot(sin, color = &output)
    wdg = Widget(["r" => r, "g" => g, "b" => b, "update" => update], output = output)
    @layout! wdg hbox(plt, vbox(:r, :g, :b, :update)) ## custom layout: by default things are stacked vertically
end


loadbutton = filepicker()
hellobutton = button("Hello!")
goodbyebutton = button("Good bye!")
ui = vbox( # put things one on top of the other
    loadbutton,
    hbox( # put things one next to the other
        pad(1em, hellobutton), # to allow some white space around the widget
        pad(1em, goodbyebutton),
    )
)
display(ui)

loadbutton = filepicker()
columnbuttons = Observable{Any}(dom"div"())

display(ui)

ui = vbox( # put things one on top of the other
    loadbutton,
    hbox( # put things one next to the other
        pad(1em, hellobutton), # to allow some white space around the widget
        pad(1em, goodbyebutton),
    )
)
display(ui)
using Colors

@manipulate for r = 0:.05:1, g = 0:.05:1, b = 0:.05:1
    HTML(string("<div style='color:#", hex(RGB(r,g,b)), "'>Color me</div>"))
end

timepicker(value::Union{Dates.Time, Observable, Nothing}=nothing)

textarea(hint=""; value="")

@manipulate throttle=.05 for λ=0:.1:5, μ=0:.1:5
    xs = range(0.0, 1.0, length = 100)
    Plots.plot(xs, x -> λ*x^2 + μ)
end

plotlyjs()
plot(x, y)

x = 1:1:10
y = 1:1:10

using Blink
using Plots
#--------------------#
a = [0 10 -4 1;1 4 -1 1;3 2 1 2;-2 -8 2 -2 ;1 -6 3 0]
m = [1;2;5;-4;1]
x = inv(a)*m


#_-----------------------#
#----------------P1----------------------#

# result is ->
#-5.0 - 1.0*a*(1 + a) + 1.0*a*(-2 + a)*(1 + a) + 4.0*(1 + a)



#------------------P2----------------------
#r1  = 1.0 + 1.0*a - 1.3333333333*a*(-2 + a)
#r2 = -0.0833333333*(-2 + a)*(1 + a) + 0.0416666667*(-3 + a)*(-2 + a)*(1 + a) + 0.3333333333*(1 + a)
#r3 =  -2.0 + 1.5*a


#-----------------P3------------------------
#r4 = 3.0 + 1.0*(1 + a)*(-1 + a) - 1.0*(1 + a)

#------------P4-----------------------------
#r5 = 8.0 - 2.0*(2 + a)

#----------P5----------------
# r6 = 1.1023e-06*(-3 + a)*(-2 + a)*(-5 + a)*(-4 + a)*(-1 + a)*(-10 + a)*(-7 + a)*(-6 + a)*(-9 + a)*(-8 + a)


#--------------------------
y = log.(x.+1)

NewtonDissection(x,y)
langrage1(x,y,0,1)

g =b-> 0.5*(b*(b+1))^(1/2)*(-0.5)/(b*(b-1))-0.75

SecantMethod(g,0,1)

x = 0:0.001:1
plot(x,g.(x),linewidth=3)
using Plots

x = [1 1.05 1.10 1.15]
y = [0.1924 0.2414 0.2933 0.3492]

f=x->log(10,tan.(x))
g = a->1.4666666667*a^3 - 4.040000000105*a^2 + 4.6383333334435*a - 1.8726000000385

f1 = derivative(f)
f2=derivative(f1)
f3=derivative(f2)
f4 = derivative(f3)
