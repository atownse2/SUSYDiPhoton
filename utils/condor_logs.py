import os

condor_dir = '/scratch365/atownse2/SUSYDiPhoton/condor'

def clear_logs(job_names=None):

    os.system(f"rm {condor_dir}/log/*" )
    os.system(f"rm {condor_dir}/out/*" )
    os.system(f"rm {condor_dir}/err/*" )


if __name__ == '__main__':
    
    import argparse
    parser = argparse.ArgumentParser()

    parser.add_argument('--clear_logs', action='store_true', help='Clear all condor logs')

    args = parser.parse_args()

    if args.clear_logs:
        clear_logs()
        print('Cleared all condor logs')