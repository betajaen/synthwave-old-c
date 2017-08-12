function tutorial(_name)
  project(_name)
    kind           "ConsoleApp"
    language       "C"
    objdir         "_build"
    flags           { "WinMain", "FatalWarnings" }
    defines         { }
    links           { "SDL2", "SDL2main", "SDL2_net" }
    includedirs     { "..", "ref/SDL2/include",  "ref/SDL2_net/include", "ref/" }
    libdirs         { "..", "ref/SDL2/lib/x86/", "ref/SDL2_net/lib/x86/",  }
    files           { "../synthwave.h", "../synthwave.c", "../examples/" .. _name .. ".c" }
end

solution "Synthwave"

  configurations      { "Debug", "Release" }

  targetdir           "bin"
  debugdir            "bin"

  configuration "Debug"
    defines           { "DEBUG",  }
    flags             { "Symbols" }

  configuration "Release"
    defines           { "NDEBUG"}
    flags             { "Optimize" }

--------------------------------------------------------------------------

group "tutorials"

--------------------------------------------------------------------------
