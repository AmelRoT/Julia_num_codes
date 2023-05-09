f(n) = n^n/(factorial(n))
println(f(5))
t = f.(1:1:17)
for i in t
    println(i)
end
println("===========================")

function f(n)
 p = 1
for i in 1:(n-1)

p = p*n/(n-i)

end
return p

end

t = f.(1:1:700)

for j in 1:length(t)
println(t[j])

end
println("=================")

t1 = 1:1:200
for i in t1
println(i)

end
