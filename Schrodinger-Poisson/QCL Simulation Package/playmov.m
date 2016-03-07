load('wfmovie'); 
ctr = 1; 
flag = 1; 
while flag
    k=waitforbuttonpress;
    [Img,MAP] = frame2im(Mov(ctr)); 
    imshow(Img);
    title(num2str(bias(ctr))); 
    if strcmp(get(gcf,'currentcharacter'),'a');
        if(ctr > 1)
            ctr = ctr -1;
        end; 
    elseif strcmp(get(gcf,'currentcharacter'),'d')
        if(ctr <length(Mov) ) 
            ctr = ctr +1;
        end; 
    elseif strcmp(get(gcf,'currentcharacter'),'x')
        flag = 0; 
    else
        display('Unrecognized key. Press a for prev frame, d for next frame and x to exit!');
    end
    
end