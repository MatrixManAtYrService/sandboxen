import conducto as co
from inspect import cleandoc

env = {"HEROKU_API_KEY": "88d1c57c-c074-4333-9004-56f1b6b32e11"}


def pipeline() -> co.Parallel:
    root = co.Parallel()
    root["one"] = co.Exec(
        cleandoc(
            """
            docker run --rm \\
                -e HEROKU_API_KEY='88d1c57c-c074-4333-9004-56f1b6b32e11' \\
                dickeyxxx/heroku-cli \\
                heroku apps
            """
        ),
        requires_docker=True,
        image="docker:latest",
    )
    root["two"] = co.Exec("heroku apps", env=env, image="dickeyxxx/heroku-cli")
    return root


if __name__ == "__main__":
    co.main(default=pipeline)
