import conducto as co


def func() -> co.Exec:
    return co.Exec("echo $FOO", env={"FOO": "foo"})


def main() -> co.Serial:
    with co.Serial() as node:

        node["one"] = func()

        # doesn't work (but could)
        # node["two"] = func().set(env={"FOO": "bar"})

        # works
        node["two"] = func()
        node["two"].set(env={"FOO": "bar"})

    return node


if __name__ == "__main__":
    co.main(default=main)
