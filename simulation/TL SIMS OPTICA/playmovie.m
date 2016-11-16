load('Movie'); 
ctr = 1; 
flag = 1; 
[sort_biases,idx] = sort(biases); 
biases = biases(idx); 
MOV = MOV(idx);

while flag
    k=waitforbuttonpress;
    [Img,MAP] = frame2im(MOV(ctr)); 
    imshow(Img);
    title(num2str(biases(ctr))); 
    if strcmp(get(gcf,'currentcharacter'),'a');
        if(ctr > 1)
            ctr = ctr -1;
        end; 
    elseif strcmp(get(gcf,'currentcharacter'),'d')
        if(ctr <length(MOV) ) 
            ctr = ctr +1;
        end; 
    elseif strcmp(get(gcf,'currentcharacter'),'x')
        flag = 0; 
    else
        display('Unrecognized key. Press a for prev frame, d for next frame and x to exit!');
    end
    
end