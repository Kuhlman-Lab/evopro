import multiprocessing as mp
import collections
from typing import Sequence, Union

from functools import partial
import numpy as np

class Distributor:
    """This class will distribute work to sub-processes where
    the same function is run repeatedly with different inputs.
    The Distributor is initialized with the number of worker
    processes and a function, f_init, whose job it is to create
    a function "f" that will be run repeatedly by the worker
    processes.
    """

   
    def __init__(self, n_workers, f_init, arg_file, lengths):
        """
        Construct a Distributor that manages n_workers sub-processes.
        The distributor will give work to its sub-processes in the
        form of function inputs.

        f_init should be a function and return a function. It will
        be called by each sub-process once as that process gets started.
        It should return the worker function that will do the heavy
        lifting.
        """
        self.n_workers = n_workers
        self.lock = mp.Lock()
        self.qs_out = [mp.Queue() for _ in range(n_workers)]
        self.q_in = mp.Queue()
        self.processes = [
            mp.Process(
                target=Distributor._worker_loop,
                args=(
                    f_init,
                    i,
                    self.lock,
                    self.qs_out[i],
                    self.q_in,
		    arg_file,
		    lengths)
            )
            for i in range(n_workers)
        ]
   
        for p in self.processes:
            p.start()


    def spin_down(self):
        """When all work is done, send out a spin-down signal to all the
        subprocesses and join them
        """
        for i in range(self.n_workers):        
            # spin down the worker
            self.qs_out[i].put((False, None))
        for i in range(self.n_workers):
            self.processes[i].join()
           

    def churn(self, work_list):
        """Process the work in the work list, farming out work
        to the subprocesses
        """
        work_queue = collections.deque(work_list)
        n_jobs = len(work_queue)
        job_ind_for_worker = [-1] * self.n_workers
        job_output = [None] * n_jobs


        count_job = 0
        count_completed = 0
       
        # step 1: put all the work we already have into the work queues
        if self.n_workers < n_jobs:
            for i in range(self.n_workers):
                w = work_queue.popleft()
                self.qs_out[i].put((True, w))
                job_ind_for_worker[i] = count_job
                count_job += 1
        else:
            count = 0
            while len(work_queue) > 0:
                w = work_queue.popleft()
                self.qs_out[count].put((True, w))
                job_ind_for_worker[count] = count_job
                count += 1
                count_job += 1

        while len(work_queue) > 0:
            proc_id, val = self.q_in.get()
            count_completed += 1
            #print("work_list loop count completed:", count_completed)
            job_ind = job_ind_for_worker[proc_id]
            assert job_ind != -1
            job_output[job_ind] = val
   
            w = work_queue.popleft()
            self.qs_out[proc_id].put((True, w))
            job_ind_for_worker[proc_id] = count_job
            count_job += 1

        while count_completed < n_jobs:
            proc_id, val = self.q_in.get()
            #print("wait-for-jobs to finish loop. count completed:", count_completed)
            count_completed += 1
            job_ind = job_ind_for_worker[proc_id]
            assert job_ind != -1
            job_ind_for_worker[proc_id] = -1
            job_output[job_ind] = val

        return job_output

           
    @staticmethod
    def _worker_loop(f_init, proc_id, lock, q_in, q_out, arg_file, lengths):

        f = f_init(proc_id, arg_file, lengths)
   
        is_job, val = q_in.get()
        while is_job:
            #print(val)
            result = f(val)
       
            lock.acquire()
            try:
                q_out.put((proc_id, result))
            finally:
                lock.release()
            is_job, val = q_in.get()
        print("spinning down worker", proc_id)


if __name__=="__main__":
    print("no main functionality")
