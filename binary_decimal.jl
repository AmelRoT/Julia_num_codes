# ----------Decimal to binary conversion (not only integers)------------#

#Enter the value of the number (decimal or integer value ) ->

decimal_number = 9418241.1298512;
display("Decimal Value is : ")
display(decimal_number)
# we may use the integrated function called trunc to convert the float or double value into int.

prior_to_comma_part = trunc(Int,decimal_number); # prior to comma(numbers)

after_the_comma_part = decimal_number-prior_to_comma_part; # after the comma(numbers)


#------------------------prior to the comma part-----------------------------

m = 10000; # the purpose is to mimic the while loop (not interested in using while loop here)
value_prior = prior_to_comma_part; #takes the value that we will manipualte and find binary form of.
var1 = 0;#need to initilize the matrix in order to store the values
var2 = 0; #initilize the stroing place for prior to decimal point numbers;
for i in 1:m

    if(mod(value_prior,2)==0 && value_prior!=0) # mod function takes the remainder of the value divided by 2
    var1 = [var1 0];
    value_prior = trunc(Int,value_prior/2)
    #display(value_prior)
    end

    if(mod(value_prior,2)==1 && value_prior!=0) # mod function takes the remainder of the value divided by 2
    var1 = [var1 1];
    value_prior = trunc(Int,value_prior/2);

    end

    if(value_prior==0)
        break;
    end

end

    var2 = var1[length(var1)];
    for i in 1:length(var1)-2

    var2 = [var2 var1[length(var1)-i]]

    end

#------------------------------end of the prior to the comma part -------------------------#

#------------------------------After the comma part --------------------------------------#

value_after = after_the_comma_part;
var3 = 0;
counter = 1 ;
for j in 1:1000
    if(value_after!=0)
        var3 =[var3 trunc(Int,value_after*2)];
        value_after = mod(value_after*2,1)
    end

    if(value_after==0 || counter == 20) # in julia 0.0000000 == 0 , that is why we don't need to use or value_after ==0.0 etc

        break;

    end
    counter = counter+1; # we will limit the after the comma part to twenty numbers.
end
if(after_the_comma_part!=0)
var4 = var3[2]
end

for h in 3:length(var3)
    var4 = [var4 var3[h]];
end

#-----------------------------Adding everything together-----------------------------------
if(after_the_comma_part!=0)
display("The binary value for the provided decimal number is : ")
display(string(var2)*"."*string(var4))
end
if(after_the_comma_part==0)
    display("The binary value for the provided decimal number is : ")
    display(string(var2))
end
# ----------Binary to decimal conversion (not only integers)------------#
# [1 1 0 1 0] ˛& [0 1 0 0 0 0 ] == 1 1 0 1 0.0 1 0 0 0 0 || we will store the binary number in array

binary_number_prior_to_comma = [1 1 0 0 0 1 1 1 1 1 0 1]; #prior to comma part, needs spacing
binary_number_after_the_comma = [0 1 1 0 0 1 1 0 1]; #after the comma part, needs spacing
display(string(binary_number_prior_to_comma)*"."*string(binary_number_after_the_comma))
value_decimal_prior = 0; #storing value
for i in 1:length(binary_number_prior_to_comma)

    value_decimal_prior = value_decimal_prior + 2^(length(binary_number_prior_to_comma)-i)*binary_number_prior_to_comma[i];
#    display(value_decimal_prior)
end

value_decimal_after = 0; #storing value
for i in 1:length(binary_number_after_the_comma)

    value_decimal_after = value_decimal_after + 2.0^(-i)*binary_number_after_the_comma[i];
    #display(value_decimal_after)
end

#The final decimal value is ->
decimal_number_value_final = value_decimal_prior+value_decimal_after;
    display("The decimal value for the provided binary number is : ")
    display(decimal_number_value_final)
