function printINFO( info )

        display([info.settings.name ' ' info.cavity ' cavity simulation in a ' num2str(info.Ltot) '(mm) long ridge']);
        display(['iteration Nr. = ' num2str(info.iter_ctr) ' @ RT = ' num2str(info.RT) '/' num2str(info.settings.simRT) ])
        display(['Grid size = ' num2str(info.N) ]);
        display(['Current max intensity: ' num2str(info.maxInt) ]);

end

