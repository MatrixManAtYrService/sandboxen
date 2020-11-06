from inspect import cleandoc
from calendar import monthrange
from datetime import datetime
from textwrap import indent, dedent, wrap
import conducto as co
import os


def pipeline() -> co.Serial:

    # defer node definition until the first node runs
    root = co.Lazy(nodes_for_this_month)

    # conducto installs the dependencies into its image
    root.image = co.Image(
        copy_url="https://github.com/MatrixManAtYrService/sandboxen",
        copy_branch="master",
        path_map={".": "./fortune_witherror"},
        reqs_py=["conducto", "sh"],
        reqs_packages=["fortune"],
    )

    return root


def nodes_for_this_month() -> co.Parallel:
    """
    Nodes defined as Lazy become two nodes in a pipeline instance:

    - Generate
    - Execute

    The Generate node calls this function.
    The Execute node runs whichever node is returned by this function.
    """

    # delete to fix
    raise Exception("Some exceptions just want to watch the world burn")

    # The image parameters that appear in `reqs_py` and `reqs_packages` are
    # depenencies of this function. But the pipeline launcher doesn't need them.
    #
    # Import them inside the function to reduce external dependency.

    os.environ["PATH"] = ":".join([os.environ["PATH"], "/usr/games"])
    from sh import fortune

    now = datetime.now()
    parent = co.Parallel()
    for i in range(monthrange(now.year, now.month)[1]):

        date = f"{now.year}-{now.month}-{i + 1}"

        fortune_str = indent(fortune().stdout.decode(), prefix=16 * " ")

        cmd = cleandoc(
            f"""
            echo "About {date} the spirits say:"
            cat << EOF

                {fortune_str[16:]}
            EOF"""
        )

        parent[date] = co.Exec(cmd)

    return parent


if __name__ == "__main__":
    co.main(default=pipeline)
