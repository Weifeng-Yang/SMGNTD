% timerun(end)>timeend : Where 'timeend' refers to the termination running time

function stop=stopcheck(Res,timerun,stopindex)
    stop=0;
     if(stopindex==1)
        if(timerun(end)>30)
        stop=1;
        end
    elseif(stopindex==2)
        if(timerun(end)>60)
        stop=1;
        end
     elseif(stopindex==4)
        if(Res<1e-5)
        stop=1;
        end
    end
end