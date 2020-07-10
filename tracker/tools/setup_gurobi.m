function setup_gurobi(configurations)
%TODO throw error message
try
    envValue = getenv('GRB_LICENSE_FILE');
    if isempty(envValue)
        disp(['SETTING GRB LICENSE FILE VARIABLE to ',configurations.Solver.gurobiLic])
        setenv('GRB_LICENSE_FILE',configurations.Solver.gurobiLic);
    end
catch EM
    warning('Error while registering Gurobi')
    rethrow(EM)
end
disp('Gurobi successfully registered')

end