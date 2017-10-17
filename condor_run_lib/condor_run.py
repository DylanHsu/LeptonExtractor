import os
import time
import shutil
import subprocess
import re
# cannot use htcondor python binding because CMS python environment does not have the libraries

libdir = os.path.dirname(os.path.realpath(__file__))

class CondorRun(object):

    def __init__(self, executable):
        self.executable = executable
        self.transfer_exec = False
        self.hold_on_fail = False
        self.requirements = ''
        self.arch = 'X86_64'
        self.group = 'group_t3mit'
        self.aux_input = []
        self.min_memory = 1500
        self.num_repeats = 1

        self.job_args = []
        self.pre_args = ''
        self.post_args = ''

        self.job_names = []

        self.logdir = '/tmp'
        self.clear_log = False

        self.last_submit = [] # [(cluster id, job name)] of last submit if job_names are set

    def submit(self, name = ''):
        if name:
            logdir = self.logdir +'/' + name
        else:
            logdir = self.logdir + '/' + str(int(time.time()))

        if self.clear_log:
            try:
                shutil.rmtree(logdir)
            except:
                pass

        try:
            os.mkdir(logdir)
        except OSError:
            pass

        if len(self.job_args):
            job_args = list(self.job_args)
        else:
            job_args = ['']

        use_job_names = (len(set(self.job_names)) == len(job_args))

        with open(logdir + '/env.sh', 'w') as envfile:
            for key, value in os.environ.items():
                if key.startswith('BASH_FUNC_'):
                    envfile.write(key.replace('BASH_FUNC_', '') + ' ' + value[2:] + '\n')
                    envfile.write('export -f ' + key.replace('BASH_FUNC_', '').replace('()', '') + '\n')
                else:
                    envfile.write('export ' + key + '="' + value.replace('"', '\\"') + '"\n')

        with open(logdir + '/jobs.dat', 'w') as jobs_file:
            jobs_file.write('[EXECUTABLE] ' + os.path.realpath(self.executable) + '\n')
            jobs_file.write('[NJOBS_PER_ARG] ' + str(self.num_repeats) + '\n')
            jobs_file.write('[ARGUMENTS]\n')
            for ijob, job_arg in enumerate(job_args):
                if use_job_names:
                    jobs_file.write('%s: %s\n' % (self.job_names[ijob], job_arg))
                else:
                    jobs_file.write('%d: %s\n' % (ijob, job_arg))

        jdl = []
        jdl.append(('executable', libdir + '/condor-run.exec'))
        jdl.append(('universe', 'vanilla'))
        jdl.append(('should_transfer_files', 'YES'))
        jdl.append(('input', '/dev/null'))
        requirements = 'Arch == "%s"' % self.arch
        if self.requirements:
            requirements += ' && (%s)' % self.requirements
        jdl.append(('requirements', requirements))
        jdl.append(('rank', '32 - TARGET.SlotID')) # evenly distribute jobs across machines

        input_files = [logdir + '/env.sh'] + self.aux_input
        if self.transfer_exec:
            input_files.append(os.path.realpath(self.executable))
        jdl.append(('transfer_input_files', ','.join(input_files)))

        jdl.append(('transfer_output_files', '""'))
        jdl.append(('accounting_group', self.group))
        jdl.append(('use_x509userproxy', 'True'))
        jdl.append(('x509userproxy', '/tmp/x509up_u{0}\n'.format(os.getuid())))

        jdl.append(('request_memory', self.min_memory))
        if self.hold_on_fail:
            jdl.append(('on_exit_hold', '(ExitBySignal == True) || (ExitCode != 0)'))
       
        if use_job_names:
            jdl.append(('log', logdir + '/$(Cluster).$(JobName).log'))
            jdl.append(('output', logdir + '/$(Cluster).$(JobName).out'))
            jdl.append(('error', logdir + '/$(Cluster).$(JobName).err'))
        else:
            jdl.append(('log', logdir + '/$(Cluster).$(Process).log'))
            jdl.append(('output', logdir + '/$(Cluster).$(Process).out'))
            jdl.append(('error', logdir + '/$(Cluster).$(Process).err'))

        if self.transfer_exec:
            arg = '"./' + os.path.basename(self.executable)
        else:
            arg = '"' + re.sub('^/export', '', os.path.realpath(self.executable))

        arg += ' ' + self.pre_args + ' $(JobArgs) ' + self.post_args
        if self.num_repeats > 1:
            arg += ' $(Step)'

        arg += '"'

        jdl.append(('arguments', arg))

        jdl_text = ''.join(['%s = %s\n' % (key, str(value)) for key, value in jdl])
        jdl_text += 'queue %d ' % self.num_repeats
        if use_job_names:
            jdl_text += 'JobName,'
        jdl_text += 'JobArgs from (\n'

        for ijob, job_arg in enumerate(job_args):
            if use_job_names:
                jdl_text += self.job_names[ijob] + ', '

            if ijob == 0 and job_arg.strip() == '':
                # condor bug? cannot handle empty argument in the first line
                jdl_text += '_DUMMY_\n'
            else:
                jdl_text += job_arg + '\n'

        jdl_text += ')\n'

        with open(logdir + '/jdl', 'w') as jdl_file:
            jdl_file.write(jdl_text)

        proc = subprocess.Popen(['condor_submit'], stdin = subprocess.PIPE, stdout = subprocess.PIPE, stderr = subprocess.STDOUT)
        out, err = proc.communicate(jdl_text)
        print out.strip()

        matches = re.search('job\(s\) submitted to cluster ([0-9]+)\.', out)
        if matches:
            cluster_id = int(matches.group(1))
            if use_job_names:
                self.last_submit = []
                for ijob, job_name in enumerate(self.job_names):
                    self.last_submit.append(('%d.%d' % (cluster_id, ijob), job_name))
        else:
            print 'No cluster id found!'
            cluster_id = 0
    
        print 'Logdir is', logdir
 
        return cluster_id
