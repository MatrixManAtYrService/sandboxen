import conducto as co
from pathlib import Path
import sys


def fail_then_pass(key: str, times: int):

    file = Path(f'/conducto/data/pipeline/{key}')

    # increment file counter
    if not file.exists():
        with open(file, 'w') as f:
            f.write('1')
        ct = 1
    else:
        with open(file, 'r') as f:
            ct = int(f.read())
        ct += 1
        with open(file, 'w') as f:
            f.write(str(ct))

    if ct > times:
        sys.exit(0)
    else:
        sys.exit(1)


def main() -> co.Serial:

    retry = co.Exec(fail_then_pass, "retry", 2)
    retry.on_error(co.callback.retry(3))

    retry_2 = co.Exec(fail_then_pass, "retry2", 3)
    retry_2.on_error(co.callback.retry(2))

    retry_then_skip = co.Exec(fail_then_pass, "retry_then_skip", 3)
    retry_then_skip.on_error(co.callback.retry_then_skip(2))

    retry_then_skip_2 = co.Exec(fail_then_pass, "retry_then_skip", 2)
    retry_then_skip_2.on_error(co.callback.retry_then_skip(3))

    retry_with_double_mem = co.Exec(fail_then_pass, "retry_with_double_mem", 2)
    retry_with_double_mem.on_error(co.callback.retry_with_double_mem(3))

    handle_memory_errors = co.Exec(fail_then_pass, "retry_with_double_mem", 0)
    handle_memory_errors.on_error(co.callback.handle_memory_errors())

    skip_some_errors = co.Serial(stop_on_error=False)
    skip_some_errors['pass'] = co.Exec('echo hi')
    skip_some_errors['fail1'] = co.Exec('echo hi | grep foo')
    skip_some_errors['fail2'] = co.Exec('echo hi | grep bar')
    skip_some_errors.on_error(co.callback.skip_some_errors(2))

    with co.Serial(image=co.Image(copy_dir=".", reqs_py=["conducto"]), stop_on_error=False) as node:
        node["retry"] = retry
        node["retry_2"] = retry_2
        node["retry_then_skip"] = retry_then_skip
        node["retry_then_skip_2"] = retry_then_skip_2
        node["retry_with_double_mem"] = retry_with_double_mem
        node["skip_some_errors"] = skip_some_errors
        node["handle_memory_errors"] = handle_memory_errors

    node.on_done(co.callback.email(to="mrixman@conducto.com"))

    return node


if __name__ == "__main__":
    co.main(default=main)
