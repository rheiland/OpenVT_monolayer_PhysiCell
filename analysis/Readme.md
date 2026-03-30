
```
(base) M1P~/git/monolayer/analysis$ g++-14 metrics.cpp 
In file included from /opt/homebrew/Cellar/gcc/14.2.0/include/c++/14/bits/stl_algobase.h:64,
                 from /opt/homebrew/Cellar/gcc/14.2.0/include/c++/14/string:51,
                 from /opt/homebrew/Cellar/gcc/14.2.0/include/c++/14/bits/locale_classes.h:40,
                 from /opt/homebrew/Cellar/gcc/14.2.0/include/c++/14/bits/ios_base.h:41,
                 from /opt/homebrew/Cellar/gcc/14.2.0/include/c++/14/ios:44,
                 from /opt/homebrew/Cellar/gcc/14.2.0/include/c++/14/ostream:40,
                 from /opt/homebrew/Cellar/gcc/14.2.0/include/c++/14/iostream:41,
                 from metrics.cpp:1:
/opt/homebrew/Cellar/gcc/14.2.0/include/c++/14/bits/stl_pair.h: In instantiation of 'constexpr std::pair<typename std::__strip_reference_wrapper<typename std::decay<_Tp>::type>::__type, typename std::__strip_reference_wrapper<typename std::decay<_Tp2>::type>::__type> std::make_pair(_T1&&, _T2&&) [with _T1 = double&; _T2 = double&; typename __strip_reference_wrapper<typename decay<_Tp>::type>::__type = double; typename decay<_Tp>::type = double; typename __strip_reference_wrapper<typename decay<_Tp2>::type>::__type = double; typename decay<_Tp2>::type = double]':
metrics.cpp:36:26:   required from here
   36 |     return std::make_pair(a[1], a[0]) < std::make_pair(b[1], b[0]);
      |            ~~~~~~~~~~~~~~^~~~~~~~~~~~
/opt/homebrew/Cellar/gcc/14.2.0/include/c++/14/bits/stl_pair.h:1132:5: note: parameter passing for argument of type std::pair<double, double>' when C++17 is enabled changed to match C++14 in GCC 10.1
 1132 |     make_pair(_T1&& __x, _T2&& __y)
      |     ^~~~~~~~~
(base) M1P~/git/monolayer/analysis$ ll metrics*
-rwxr-xr-x  1 heiland  staff  345960 Oct 13 16:13 metrics*
-rw-r--r--  1 heiland  staff    5576 Oct 13 16:13 metrics.cpp
(base) M1P~/git/monolayer/analysis$ ll a.out 
-rwxr-xr-x  1 heiland  staff  346504 Oct 13 16:20 a.out*
(base) M1P~/git/monolayer/analysis$ 
```
