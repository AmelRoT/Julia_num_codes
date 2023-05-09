
a = [0 30 60 90 120 150 180 210 240 270 300 330]

b = sind.(a)

b_norm = b.+1

perc_b = (b_norm/2)*100



# with 10 deg
a = [0]
for i = 1: 35

a = [a a[i]+10]


end

b = sind.(a)
b_norm = b.+1

perc_b = (b_norm/2)*100

# with 5 deg
a = [0]
for i = 1: 71

a = [a a[i]+5]


end

b = sind.(a)
b_norm = b.+1

perc_b = (b_norm/2)*100
