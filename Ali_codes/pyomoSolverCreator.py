from pyomo.environ import *
from pyomo.opt import *
from userError import userError
import sys
sys.path.append('/Users/elieeshoa/Dropbox/Elie_Eshoa/Ali_codes/GAMETES/suppl_material/')
from globalVariables import *

def pyomoSolverCreator(optSolverName, max_threads_num = 1, **solver_options):
    """
    Creates a pyomo solver object and assigns the solver options

    INPUTS:
    ------
     optSolverName: A string containing the solver name
    solver_options: Additional solver options. These are to a dictionary whose keys are 
                    are the names of the arguments and values are the values 

    OUTPUTS:
    -------
    pymoSolverObject: Pyomo solver object with solver options assigned
    """    
    if not isinstance(max_threads_num,int):
        raise TypeError('An integer expected for max_threads_num. An object of type {} was provided instead'.format(type(max_threads_num)))
    elif max_threads_num < 1 or max_threads_num > 32:
        raise ValueError('An integer between 1 and 32 is expected for max_threads_num')

    # Pyomo solver object
    pymoSolverObject = SolverFactory(optSolverName)

    # - Set some default solver options -        
    if optSolverName.lower() == 'cplex':
        # Memory
        # https://www.ibm.com/support/knowledgecenter/en/SS9UKU_12.4.0/com.ibm.cplex.zos.help/Parameters/topics/WorkMem.html
        pymoSolverObject.options["workmem"]=2500

        # Feasbility tolerance (eprhs). Defaul = 1e-6
        pymoSolverObject.options["simplex_tolerances_feasibility"]=1e-9

        # Optimality tolerance (epopt). Default = 1e-6
        pymoSolverObject.options["simplex_tolerances_optimality"]=1e-9

        # Relative MIP optimality gap, Default = 1e-4
        pymoSolverObject.options["mip_tolerances_mipgap"] = 1e-9

        # Absolute MIP optimality gap, Default = 1e-6
        pymoSolverObject.options["mip_tolerances_absmipgap"] = 1e-10

        # Integrality tolerance (epint). Default = 1e-5
        pymoSolverObject.options["mip_tolerances_integrality"] = mip_integrality_tol

        # MIP strategy variable select (varsel). Default = 0
        # https://www.ibm.com/support/knowledgecenter/en/SS9UKU_12.4.0/com.ibm.cplex.zos.help/Parameters/topics/VarSel.html
        # NOTE: Sometimes using a value of 3 can significantly increase the runtime (my personal experience) 
        #pymoSolverObject.options["mip_strategy_variableselect"]=3

        # Bound strengthening indicator (bndstrenind). Default = -1 (let cplex decides) 
        # https://www.ibm.com/support/knowledgecenter/ru/SSSA5P_12.6.1/ilog.odms.cplex.help/CPLEX/Parameters/topics/BndStrenInd.html
        pymoSolverObject.options["preprocessing_boundstrength"]=1

        # Number of threads (core) to use (this shouldn't be more than the value specified in your job
        # https://www.ibm.com/support/knowledgecenter/en/SS9UKU_12.5.0/com.ibm.cplex.zos.help/Parameters/topics/Threads.html
        # by #$-pe omp N otherwise it may use more cores than N and the job gets killed
        # by the system as it uses more cores than requested 
        pymoSolverObject.options["threads"] = max_threads_num

        # Directory for temporary fiiles
        #pymoSolverObject.options["workdir"] = pyomo_tmp_dir
      
    elif optSolverName.lower() in 'gurobi':
        # Memory (in Gb). Default: Infinity
        # http://www.gurobi.com/documentation/7.0/refman/nodefilestart.html
        pymoSolverObject.options["NodefileStart"]=2.5

        # Feasbility tolerance. Defaul = 1e-6
        pymoSolverObject.options["FeasibilityTol"]=1e-9

        # Optimality tolerance. Defaul = 1e-6
        pymoSolverObject.options["OptimalityTol"]=1e-9

        # Relative MIP optimality gap, Default = 1e-4
        pymoSolverObject.options["MIPGap"]=1e-9

        # Absolute MIP optimality gap, Default = 1e-10
        pymoSolverObject.options["MIPGapAbs"]=1e-10

        # Integrality tolerance (epint). Default = 1e-5
        pymoSolverObject.options["IntFeasTol"]= mip_integrality_tol

        # Branch variable selection strategy . Default = -1 (automatic)
        # http://www.gurobi.com/documentation/7.5/refman/varbranch.html
        # NOTE: Sometimes using a value of 3 can significantly increase the runtime 
        #pymoSolverObject.options["VarBranch"]=3

        # LP method used to solve sifting sub-problems. Default = -1 (automatic)
        # Setting this parameter to 2 (barrier method) somethimes helps when the solver terminates 
        # with a "sub-optimal solution"
        # http://www.gurobi.com/documentation/7.5/refman/siftmethod.html#parameter:SiftMethod 
        pymoSolverObject.options["SiftMethod"]=2

        # Number of cores to use for parallel computations (defaul 0 = all cores) 
        # NOTE: If you use the default value of 0 and submit the code on clusters, it uses all available cores 
        # (e.g., 16) and the job is usually killed unless you request the maximum number of processors in your code
        # (e.g., #$-pe omp 16). So, it's best to set this to a value that is reasonable (in terms of your job's
        # wait time in the queue based on the number of cores you are requesting).
        # The value of this parameter shoudl always be less than that used in #$-pe omp 
        pymoSolverObject.options["Threads"] = max_threads_num 

    elif optSolverName.lower() in 'gurobi_ampl':
        # Memory (in Gb). Default: Infinity
        pymoSolverObject.options["NodefileStart"]=2.5

        # Feasbility tolerance. Defaul = 1e-6
        pymoSolverObject.options["feastol"]=1e-9

        # Optimality tolerance. Defaul = 1e-6
        pymoSolverObject.options["opttol"]=1e-9

        # Integrality tolerance (epint). Default = 1e-5
        pymoSolverObject.options["IntFeasTol"]= mip_integrality_tol

        # Branch variable selection strategy . Default = -1 (automatic)
        # NOTE: Sometimes using a value of 3 can significantly increase the runtime (my personal experience) 
        #pymoSolverObject.options["VarBranch"]=3

        # Number of cores to use for parallel computations (defaul 0 = all cores) 
        # NOTE: If you use the default value of 0 and submit the code on clusters, it uses all available cores 
        # (e.g., 16) and the job is usually killed unless you request the maximum number of processors in your code
        # (e.g., #$-pe omp 16). So, it's best to set this to a value that is reasonable (in terms of your job's
        # wait time in the queue based on the number of cores you are requesting).
        # The value of this parameter shoudl always be less than that used in #$-pe omp 
        pymoSolverObject.options["Threads"] = max_threads_num 

    else:
        raise userError('**Error! Invalid optimization solver name ...')

    #--- Set provided solver options by the user (these overwrfite tho defaults)
    for option_name in solver_options.keys():
        exec('pymoSolverObject.options["' + option_name+ '"] = ' + str(solver_options[option_name]))


    return pymoSolverObject 
