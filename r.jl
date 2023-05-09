    a = 1
for i = 2:36

a = [a i]
end
a = [a 0]

b = 1
for j = 4:3:36

b = [b a[j]]

end


c = 2
for j = 5:3:36

c = [c a[j]]

end

d = 3
for j = 6:3:36

d = [d a[j]]

end

for j = 1:5
    mon = 140
    bet1 = 50
    bet2 = 50
    bet3 = 40
    cumlSum = mon;
    mon1 = 0;
    for ii = 1:3

    val = rand(a)

    display("Value is : ")
    display(val)
    display("--------------------------------")

    if(  val == b[3] || val == b[4] || val == b[5] || val == b[6] || val == b[7] ||
         val == b[8] || val== b[9] || val == b[10] || val == b[11] || val == b[12])

         bet1 =bet1*3
         mon1 = bet1-bet2-bet3-50
         cumlSum = cumlSum+mon1

     end

    if(val == c[1] || val == c[2] || val == c[3] || val == c[4] || val == c[5] ||
        val == c[6] || val == c[7] ||
         val == c[8] || val == c[9] || val== c[10] || val == c[11] || val == c[12])

         bet2 =bet2*3
         mon1 = bet2-bet1-bet3-50
         cumlSum = cumlSum+mon1

     end

     if(val == d[1] || val == d[2] || val == d[3] || val == d[4] || val== d[5] ||
         val== d[6] || val== d[7] ||
          val == d[8] || val == d[9] || val == d[10] || val==d[11] || val==d[12])

          bet3 =(bet3/10)*36
          mon1 = bet3-bet1-bet2-40
          cumlSum = cumlSum+mon1

      end

      if(val == 0 || val == b[1] || val == b[2])
          display("LOST MONEY")
          cumlSum = cumlSum-mon
      end
    bet1 = 50
    bet2 = 50
    bet3 = 40
    display(mon1)


    end

    final = cumlSum
    display(final)
    display("------------------------------------------------------")
    display("------------------------------------------------------")
    display("------------------------------------------------------")
    display("------------------------------------------------------")

end
