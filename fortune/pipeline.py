from inspect import cleandoc
from calendar import monthrange
from datetime import datetime
from textwrap import indent
import conducto as co
import os


def nodes_for_this_month() -> co.Parallel:
    """
    This function runs in the container for the generate step.
    It returns a node to be executed as part of the execute step.
    """

    # linux utility: fortune
    # python library: sh
    # The above are not dependencies for launching this pipeline, but they
    # must be installed in the image to be referenced by this function.

    os.environ['PATH'] = ':'.join([os.environ['PATH'], "/usr/games"])
    from sh import fortune

    now = datetime.now()
    parent = co.Parallel()
    for i in range(monthrange(now.year, now.month)[1]):

        date = f"{now.year}-{now.month}-{i + 1}"

        fortune_str = fortune()

        cmd = cleandoc(
            f"""
            echo "About {date} the spirits say:"
            cat << EOF
            {indent(fortune_str, prefix='            ')}
            EOF"""
        )

        parent[date] = co.Exec(cmd)

    return parent

# copy_dir places this file in the image so that
# the above function can be found when the Lazy node runs
img = co.Image(copy_dir=".",
               reqs_py=["conducto", "sh"],
               reqs_packages=["fortune"])


def make_pipeline() -> co.Serial:
    root = co.Serial(image=img)
    root['fortune'] = co.Lazy(nodes_for_this_month)
    return root


if __name__ == "__main__":
    co.Image.share_directory("fortune", ".")
    co.main(default=make_pipeline)
