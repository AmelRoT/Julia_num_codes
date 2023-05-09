println("Enter your number : \n")
num = readline();
num = parse(Int64, num);
println("The number is equal to :  ", num)
i = 0
l = [10, 20, 30, 40, 50, 60]
for i in l
    println(1)
end

a = 1;
m = 1;
for i = 1:5
    for j = 1:m
        print(a, " ")
    end
    println()
    a = a + 1
    m = m + 1
end

i = 0;
while (i < 5)
    if (i == 4)
        break
    end
    i = i + 1
    println(i)
end

x = 3;
x = -10
y = x^2 - 3 * x - 4;
y1 = 2 * x - 3;
a = y / y1
println(a)

for i = 1:20
    x1 = x - (y / y1)
    println(x1)
    x = x1
    y = x^2 - 3 * x - 4
    y1 = 2 * x - 3
end

x = collect(-10:0.001:10)
y = x .^ 2 - 3 * x .- 4;
using Plots
plotly()
plot(x, y, linewidth = 4)
plot!([4.0], [0.0], markershape = :circle, color = "red", markersize = 10)
plot!([-1.0], [0.0], markershape = :circle, color = "red", markersize = 10)
