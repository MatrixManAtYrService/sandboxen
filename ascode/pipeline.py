import conducto as co

IMG = co.Image("python:3.9", copy_dir=".", reqs_docker=True)


def main() -> co.Serial:
    with co.Serial(image=IMG, requires_docker=True) as root:
        with co.Parallel(name="Init") as init:
            init["Build"] = co.Exec("docker build .")
            init["Lint"] = co.Exec("black --check .")
            init["Unit Test"] = co.Exec("python unit_test.py")
        root["Deploy"] = co.Exec("bash deploy_aws.sh")
        root["Integration Test"] = co.Exec("bash int_test.sh")
    return root


if __name__ == "__main__":
    co.main(default=main)
