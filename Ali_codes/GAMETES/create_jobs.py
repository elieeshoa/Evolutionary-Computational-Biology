job = \
"""#!/bin/bash

# Job parameters
#BSUB -J job_{SIZE}
#BSUB -o /PHShome/ee869/Desktop/Elie-Zomorrodi/Ali_codes/GAMETES/jobs/output/job_{SIZE}-%J.out
#BSUB -N

# The following line says to write the std output into the file to see the 
# progress of the code in <desired_job_name>
LSB_STDOUT_DIRECT=y

# The following line specifies a specific queue:
# (Alternatively you can use: # bsub -q normal < my-script.lsf)
#BSUB -q normal

# --- User commands ----
# Activate conda (ignore if you don?~@~Yt have anaconda installed)
. /apps/software-2.12/Anaconda3/5.2.0/etc/profile.d/conda.sh

# Load conda (ignore if you don?~@~Yt have anaconda installed)
conda activate /PHShome/ee869/.conda/envs/manuscript

# Start time
start_time=$(date)
echo -e "\\n*** job_{SIZE} started at " $start_time "\\n"

#*** Put your commands here ***
python /PHShome/ee869/Desktop/Elie-Zomorrodi/Ali_codes/GAMETES/run_games.py --size {SIZE}
"""

# write the above string to a file called job_{SIZE}.lsf into teh directory
# jobs/scripts
# Note that the string {SIZE} will be replaced by the actual size of the game
# when the file is written

sizes = [20, 315, 385, 445, 500, 550, 590, 630, 670]
# loop over the sizes in reverse order
for size in sizes[::-1]:
    with open('jobs/scripts/job_{}.lsf'.format(size), 'w') as f:
        f.write(job.format(SIZE=size))
    print('bsub < jobs/scripts/job_{}.lsf'.format(size))


job_memory = \
"""#!/bin/bash

# Job parameters
#BSUB -J job_{SIZE}
#BSUB -o /PHShome/ee869/Desktop/Elie-Zomorrodi/Ali_codes/GAMETES/jobs/output/job_{SIZE}-%J.out
#BSUB -N

# The following line says to write the std output into the file to see the 
# progress of the code in <desired_job_name>
LSB_STDOUT_DIRECT=y

# The following line specifies a specific queue:
# (Alternatively you can use: # bsub -q normal < my-script.lsf)
#BSUB -q bigmem

# --- User commands ----
# Activate conda (ignore if you don?~@~Yt have anaconda installed)
. /apps/software-2.12/Anaconda3/5.2.0/etc/profile.d/conda.sh

# Load conda (ignore if you don?~@~Yt have anaconda installed)
conda activate /PHShome/ee869/.conda/envs/manuscript

# Start time
start_time=$(date)
echo -e "\\n*** job_{SIZE} started at " $start_time "\\n"

#*** Put your commands here ***
python /PHShome/ee869/Desktop/Elie-Zomorrodi/Ali_codes/GAMETES/run_games.py --size {SIZE}
"""

sizes = [705, 740, 775, 805, 835, 865, 895, 920, 950, 975, 1000]
for size in sizes[::-1]:
    with open('jobs/scripts/job_{}.lsf'.format(size), 'w') as f:
        f.write(job_memory.format(SIZE=size))
    print('bsub < jobs/scripts/job_{}.lsf'.format(size))



# 54373   ee869   RUN   interactiv eris2n4     dn021       /bin/bash  Jun 27 21:05
# 67461   ee869   RUN   normal     eris2n4     hn005       job_950    Jun 29 16:19
# 67462   ee869   RUN   normal     eris2n4     hn005       job_920    Jun 29 16:19
# 67463   ee869   PEND  normal     eris2n4                 job_895    Jun 29 16:19
# 67464   ee869   PEND  normal     eris2n4                 job_865    Jun 29 16:19
# 67465   ee869   PEND  normal     eris2n4                 job_835    Jun 29 16:19
# 67466   ee869   PEND  normal     eris2n4                 job_805    Jun 29 16:19
# 67467   ee869   PEND  normal     eris2n4                 job_775    Jun 29 16:19
# 67468   ee869   PEND  normal     eris2n4                 job_740    Jun 29 16:19
# 67469   ee869   PEND  normal     eris2n4                 job_705    Jun 29 16:19
# 67470   ee869   PEND  normal     eris2n4                 job_670    Jun 29 16:19
# 67471   ee869   PEND  normal     eris2n4                 job_630    Jun 29 16:19
# 67472   ee869   PEND  normal     eris2n4                 job_590    Jun 29 16:19
# 67473   ee869   PEND  normal     eris2n4                 job_550    Jun 29 16:19
# 67474   ee869   PEND  normal     eris2n4                 job_500    Jun 29 16:19
# 67475   ee869   PEND  normal     eris2n4                 job_445    Jun 29 16:19
# 67476   ee869   PEND  normal     eris2n4                 job_385    Jun 29 16:19
# 67477   ee869   PEND  normal     eris2n4                 job_315    Jun 29 16:19
# 67480   ee869   PEND  normal     eris2n4                 job_975    Jun 29 16:19
# 67481   ee869   PEND  normal     eris2n4                 job_950    Jun 29 16:19
# 67482   ee869   PEND  normal     eris2n4                 job_920    Jun 29 16:19
# 67483   ee869   PEND  normal     eris2n4                 job_895    Jun 29 16:19
# 67484   ee869   PEND  normal     eris2n4                 job_865    Jun 29 16:19
# 67485   ee869   PEND  normal     eris2n4                 job_835    Jun 29 16:19
# 67486   ee869   PEND  normal     eris2n4                 job_805    Jun 29 16:19

# bkill 67461
# bkill 67462
# bkill 67463
# bkill 67464
# bkill 67465
# bkill 67466
# bkill 67467
# bkill 67468
# bkill 67469
# bkill 67470
# bkill 67471
# bkill 67472
# bkill 67473
# bkill 67474
# bkill 67475
# bkill 67476
# bkill 67477
# bkill 67480
# bkill 67481
# bkill 67482
# bkill 67483
# bkill 67484
# bkill 67485
# bkill 67486

# 431824  ee869   RUN   normal     eris2n4     hn003       job_550    Nov  7 20:50
# 431827  ee869   RUN   normal     eris2n4     hn003       job_385    Nov  7 20:50
# 431821  ee869   RUN   normal     eris2n4     hn006       job_670    Nov  7 20:50
# 431823  ee869   RUN   normal     eris2n4     hn006       job_590    Nov  7 20:50
# 431825  ee869   RUN   normal     eris2n4     hn006       job_500    Nov  7 20:50
# 431826  ee869   RUN   normal     eris2n4     hn006       job_445    Nov  7 20:50
# 431822  ee869   RUN   normal     eris2n4     hn008       job_630    Nov  7 20:50
# 431858  ee869   RUN   bigmem     eris2n4     sna002a     job_975    Nov  7 21:26
# 431860  ee869   RUN   bigmem     eris2n4     sna002a     job_920    Nov  7 21:26
# 431861  ee869   RUN   bigmem     eris2n4     sna002a     job_895    Nov  7 21:26

# bkill 431824
# bkill 431827
# bkill 431821
# bkill 431823
# bkill 431825
# bkill 431826
# bkill 431822
# bkill 431858
# bkill 431860
# bkill 431861