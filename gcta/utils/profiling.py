# gcta_python/gcta/utils/profiling.py
import time
import psutil
import os

def get_memory_usage_mb() -> float:
    """
    Get current process memory usage in MB
    
    Returns:
    --------
    float
        Memory usage in MB
    """
    process = psutil.Process(os.getpid())
    return process.memory_info().rss / (1024 * 1024)

def get_system_memory_info():
    """
    Get system memory information
    
    Returns:
    --------
    dict
        System memory info in MB and percentage
    """
    mem = psutil.virtual_memory()
    return {
        'total_gb': mem.total / (1024**3),
        'available_gb': mem.available / (1024**3),
        'used_gb': mem.used / (1024**3),
        'usage_percent': mem.percent
    }


def profile_step(step_name, step_description):
    def decorator(func):
        def wrapper(self, *args, **kwargs):
            print(f"> {step_description}")
            
            # profiling 
            # ================================
            start_time = time.time()
            start_memory = get_memory_usage_mb()
            start_system = get_system_memory_info()
            
            result = func(self, *args, **kwargs)
            
            elapsed_time = time.time() - start_time
            end_memory = get_memory_usage_mb()
            end_system = get_system_memory_info()
            memory_delta = end_memory - start_memory
            # ================================
            
            # return profile
            # if isinstance(result, dict):
            self.profile[step_name] = {
                'step_description': step_description,
                'time': elapsed_time,
                'memory_before': start_memory,
                'memory_after': end_memory,
                'memory_delta': memory_delta,
                'system_before': start_system,
                'system_after': end_system,
            }
            
            # print profile
            if self.debug:
                print(f"""
-   time: {elapsed_time:.2f}sec 
-   memory: {memory_delta:.1f}MB 
-   System: {end_system['used_gb']:.1f}GB / {end_system['total_gb']:.1f}GB ({end_system['usage_percent']:.1f}%)
-   Available: {end_system['available_gb']:.1f}GB
                    """)
            
            return result
        return wrapper
    return decorator