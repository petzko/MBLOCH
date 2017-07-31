function printINFO( info )

        display([info.simname ' ' info.cavity ' cavity simulation in a ' num2str(info.Ltot) '(mm) long ridge']);
        display(['Pump := ' num2str(info.p)])
        display(['threshold inv: ' num2str(info.d_th)])
        display(['iteration Nr. = ' num2str(info.iter_ctr) ' @ RT = ' num2str(info.RT) '/' num2str(info.simRT) ])
        display(['Grid size = ' num2str(info.N) ]);
        display(['Current max intensity: ' num2str(info.maxInt) ]);

end

